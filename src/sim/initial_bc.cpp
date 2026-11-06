#include "sim/initial_bc.hpp"
#include "core/cpu_device.hpp"

#include <chrono>
#include <cmath>
#include <vector>

#include "core/mesh.hpp"
#include "core/params.hpp"
#include "fft/transforms.hpp"
#include "sim/physical_fields.hpp"

namespace initial_bc {

double centered_uniform(DisturbanceEngine &engine) {
    // An explicit 53-bit mapping keeps seeded tests stable across standard libraries.
    return static_cast<double>(engine() >> 11) * 0x1.0p-53 - 0.5;
}

void disturb_velocity(VectorField &velocity, double amplitude, const Params &p, const Mesh &m,
                      Transforms &transforms, DisturbanceEngine &generator) {
    if (amplitude == 0.0)
        return;
    const std::size_t count = static_cast<std::size_t>(m.NZ) * (m.NYP + 1) * m.NX;
    std::array<std::vector<double>, 3> physical;
    ScalarField spectral = make_scalar_field(m);
    for (int component = 0; component < 3; ++component) {
        physical[component].resize(count);
        cpu_device().for_each(-m.NZS, m.NZS + 1, [&](int iz) {
            for (int iy = 0; iy <= m.NYS; ++iy)
                for (int ix = 0; ix <= m.NXSP; ++ix)
                    spectral.at_z(iz, iy, ix, 0) = velocity.at_z(iz, iy, ix, component);
        });
        transforms.channel_xyzfft('B', physical[component].data(), spectral.data());
    }
    const double uc = p.parm.INDEX_GRAD_PRE == 1 ? p.parm.Re / 2.0 : 1.5;
    for (int iz = 0; iz < m.NZ; ++iz)
        for (int iy = 0; iy <= m.NYP; ++iy) {
            const double scale = uc * (1.0 - m.YC[iy] * m.YC[iy]) * amplitude;
            for (int ix = 0; ix < m.NX; ++ix) {
                const auto index = physical_index(iz, iy, ix, m.NYP + 1, m.NX);
                physical[0][index] += scale * centered_uniform(generator);
                physical[1][index] += scale * centered_uniform(generator);
                physical[2][index] += scale * centered_uniform(generator) * p.mesh.INDEX_disturb_Z;
            }
        }
    for (int component = 0; component < 3; ++component) {
        transforms.channel_xyzfft('F', physical[component].data(), spectral.data());
        cpu_device().for_each(-m.NZS, m.NZS + 1, [&](int iz) {
            for (int iy = 0; iy <= m.NYS; ++iy)
                for (int ix = 0; ix <= m.NXSP; ++ix)
                    velocity.at_z(iz, iy, ix, component) = spectral.at_z(iz, iy, ix, 0);
        });
    }
}

void initialize_velocity(VectorField &velocity, const Params &p, const Mesh &m,
                         Transforms &transforms) {
    DisturbanceEngine engine(
        static_cast<DisturbanceEngine::result_type>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()) +
        static_cast<DisturbanceEngine::result_type>(1000 * m.my_id));
    initialize_velocity(velocity, p, m, transforms, engine);
}

void initialize_velocity(VectorField &velocity, const Params &p, const Mesh &m,
                         Transforms &transforms, DisturbanceEngine &engine) {
    const std::size_t count = static_cast<std::size_t>(m.NZ) * (m.NYP + 1) * m.NX;
    std::vector<double> physical(count);
    ScalarField spectral = make_scalar_field(m);
    const double uc = p.parm.INDEX_GRAD_PRE == 1 ? p.parm.Re / 2.0 : 1.5;
    const std::array<double, 3> scale{uc * std::cos(p.THETA_RAD), p.ctrl.V_Microblow0,
                                      uc * std::sin(p.THETA_RAD)};

    velocity.fill({});
    for (int component = 0; component < 3; ++component) {
        cpu_device().for_each(0, m.NZ, [&](int iz) {
            for (int iy = 0; iy <= m.NYP; ++iy)
                for (int ix = 0; ix < m.NX; ++ix) {
                    const double profile = component == 1 ? 1.0 : 1.0 - m.YC[iy] * m.YC[iy];
                    physical[physical_index(iz, iy, ix, m.NYP + 1, m.NX)] =
                        scale[component] * profile;
                }
        });
        transforms.channel_xyzfft('F', physical.data(), spectral.data());
        cpu_device().for_each(-m.NZS, m.NZS + 1, [&](int iz) {
            for (int iy = 0; iy <= m.NYS; ++iy)
                for (int ix = 0; ix <= m.NXSP; ++ix)
                    velocity.at_z(iz, iy, ix, component) = spectral.at_z(iz, iy, ix, 0);
        });
    }
    disturb_velocity(velocity, p.ctrl.EPS_init, p, m, transforms, engine);
}

void velocity_boundary(VectorBoundary &boundary, double, const Params &p, const Mesh &m) {
    boundary.fill({});
    if (m.owns_zero_x_mode()) {
        boundary.at_z(0, 0, 0, 1) = p.ctrl.V_Microblow0;
        boundary.at_z(0, 1, 0, 1) = p.ctrl.V_Microblow0;
    }
}

} // namespace initial_bc
