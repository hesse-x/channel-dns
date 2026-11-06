#include <gtest/gtest.h>
#include <cmath>
#include <limits>

#include "sim/diagnostics.hpp"
#include "io/field_io.hpp"
#include "fft/transforms.hpp"
#include "stage3_fixture.hpp"
#include "test_support.hpp"

TEST(Diagnostics, AnalyticPlaneFluctuationsAndVolumeDivergence) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    Statistics statistics(m);
    Diagnostics diagnostics(m, t, nullptr);
    const auto n = static_cast<std::size_t>(m.NXTNZ) * (m.NYP+1);
    PhysicalVectorField velocity;
    for (auto& v : velocity) v.resize(n);
    PhysicalScalarField pressure(n);
    const double kx = 2*std::acos(-1.0)/p.parm.XL;
    const double kz = 2*std::acos(-1.0)/p.parm.ZL;
    for (int z = 0; z < m.NZ; ++z)
        for (int y = 0; y <= m.NYP; ++y)
            for (int x = 0; x < m.NX; ++x) {
                const auto i = physical_index(z,y,x,m.NYP+1,m.NX);
                // y-dependent fluctuation checks the quadrature, not only constant weights.
                velocity[0][i] = m.YC[y]*std::cos(kx*m.XC[x]);
                velocity[1][i] = std::sin(kx*m.XC[x]);
                velocity[2][i] = 2*std::cos(kz*m.ZC[z]);
                pressure[i] = 5 + .5*(velocity[0][i]*velocity[0][i] +
                    velocity[1][i]*velocity[1][i]+velocity[2][i]*velocity[2][i]);
            }
    auto spectral = make_vector_field(m);
    auto scalar = make_scalar_field(m);
    for (int c = 0; c < 3; ++c) {
        t.channel_xyzfft('F', velocity[c].data(), scalar.data());
        for (int z = -m.NZS; z <= m.NZS; ++z)
            for (int y = 0; y <= m.NYS; ++y)
                for (int x = 0; x <= m.NXSP; ++x)
                    spectral.at_z(z,y,x,c) = scalar.at_z(z,y,x,0);
    }
    const auto& planes = statistics.reduce(velocity, pressure, true, true);
    const auto value = diagnostics.evaluate(planes, spectral, true);
    EXPECT_TRUE(value.finite);
    EXPECT_NEAR(value.ekt, .5*n, 1e-8*n);
    EXPECT_NEAR(value.tke, 4.0/3.0, 1e-12); // .5*(.5/3 + .5 + 2).
    EXPECT_NEAR(value.divergence_l2, std::sqrt(kx*kx/6+2*kz*kz), 1e-11);
    // Real pressure is constant after the rotational-pressure correction.
    for (const auto& plane : planes) EXPECT_NEAR(plane.sum[3]/m.NXTNZ, 5, 1e-13);
    EXPECT_FALSE(diagnostics.evaluate(planes, spectral, false).finite);
    pressure[1] = std::numeric_limits<double>::quiet_NaN();
    const auto& invalid = statistics.reduce(velocity, pressure, false, true);
    EXPECT_FALSE(diagnostics.evaluate(invalid, spectral, true).finite);
}

TEST(Diagnostics, SpectralCopyCheckCatchesImaginaryPartDiscardedByRealTransform) {
    test_support::TempDirectory directory;
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    FieldIO io(m,t,directory.path);
    auto velocity = make_vector_field(m);
    auto pressure = make_scalar_field(m);
    pressure.at_z(0,0,0,0) = {0, std::numeric_limits<double>::infinity()};
    bool finite = true;
    io.to_physical(velocity, pressure, &finite);
    EXPECT_FALSE(finite);
}
