#include "sim/statistics.hpp"
#include "core/rank_filename.hpp"
#include "core/cpu_device.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <system_error>

#include "core/mesh.hpp"

namespace {

double update(double old, double sample, std::int64_t count) {
    return (old * static_cast<double>(count) + sample) / static_cast<double>(count + 1);
}

void ensure_open(const std::ofstream &out, const std::filesystem::path &path) {
    if (!out)
        throw std::runtime_error("cannot write statistics file: " + path.string());
}

template <typename Writer> void atomic_text(const std::filesystem::path &path, Writer writer) {
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    try {
        std::ofstream out(temporary, std::ios::trunc);
        ensure_open(out, temporary);
        writer(out);
        out.flush();
        if (!out)
            throw std::runtime_error("failed writing statistics file: " + path.string());
        out.close();
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

} // namespace

Statistics::Statistics(const Mesh &mesh)
    : m_(mesh), velocity_mean_(3 * (mesh.NYP + 1)), pressure_mean_(mesh.NYP + 1),
      velocity_second_(6 * (mesh.NYP + 1)), pressure_velocity_(3 * (mesh.NYP + 1)),
      pressure_second_(mesh.NYP + 1), velocity_third_(3 * (mesh.NYP + 1)),
      velocity_fourth_(3 * (mesh.NYP + 1)),
      plane_scratch_(mesh.NYP + 1) {}

namespace {
// Four independent accumulation lanes expose SIMD without architecture-specific intrinsics.
// Fixed x-contiguous loads and explicit moments avoid component loops and a pressure pass.
template<bool Higher, bool Validate>
PlaneMoments plane_moments(const Mesh& m, const PhysicalVectorField& velocity,
                          const PhysicalScalarField& pressure, int iy) {
    std::array<std::array<double, 4>, 20> sums{};
    bool finite = true;
    for (int iz = 0; iz < m.NZ; ++iz) {
        const auto base = physical_index(iz, iy, 0, m.NYP + 1, m.NX);
        for (int ix = 0; ix < m.NX; ix += 4)
            for (int lane = 0; lane < 4 && ix + lane < m.NX; ++lane) {
                const auto i = base + ix + lane;
                const double u = velocity[0][i], v = velocity[1][i], w = velocity[2][i];
                const double uu = u*u, vv = v*v, ww = w*w;
                const double raw_p = pressure[i];
                const double p = raw_p - 0.5*(uu + vv + ww);
                if constexpr (Validate)
                    finite &= std::isfinite(u) && std::isfinite(v) && std::isfinite(w) &&
                              std::isfinite(raw_p) && std::isfinite(p);
                sums[0][lane] += u; sums[1][lane] += v; sums[2][lane] += w;
                sums[3][lane] += p;
                sums[4][lane] += uu; sums[5][lane] += vv; sums[6][lane] += ww;
                if constexpr (Higher) {
                    sums[7][lane] += u*v; sums[8][lane] += v*w; sums[9][lane] += w*u;
                    sums[10][lane] += p*u; sums[11][lane] += p*v; sums[12][lane] += p*w;
                    sums[13][lane] += p*p;
                    sums[14][lane] += uu*u; sums[15][lane] += vv*v; sums[16][lane] += ww*w;
                    sums[17][lane] += uu*uu; sums[18][lane] += vv*vv; sums[19][lane] += ww*ww;
                }
            }
    }
    PlaneMoments result;
    result.finite = finite;
    for (std::size_t j = 0; j < sums.size(); ++j) {
        result.sum[j] = (sums[j][0] + sums[j][1]) + (sums[j][2] + sums[j][3]);
        if constexpr (Validate) result.finite &= std::isfinite(result.sum[j]);
    }
    return result;
}
} // namespace

const std::vector<PlaneMoments>& Statistics::reduce(const PhysicalVectorField& velocity,
                                      const PhysicalScalarField& pressure,
                                      bool higher_moments, bool validate) {
    const std::size_t expected = static_cast<std::size_t>(m_.NZ) * (m_.NYP + 1) * m_.NX;
    if (pressure.size() != expected || velocity[0].size() != expected ||
        velocity[1].size() != expected || velocity[2].size() != expected)
        throw std::runtime_error("physical field size mismatch in statistics");
    cpu_device().for_each(0, m_.NYP + 1, [&](int iy) {
        if (higher_moments) {
            plane_scratch_[iy] = validate ? plane_moments<true, true>(m_, velocity, pressure, iy)
                                         : plane_moments<true, false>(m_, velocity, pressure, iy);
        } else {
            plane_scratch_[iy] = validate ? plane_moments<false, true>(m_, velocity, pressure, iy)
                                         : plane_moments<false, false>(m_, velocity, pressure, iy);
        }
    });
    return plane_scratch_;
}

void Statistics::accumulate(const PhysicalVectorField& velocity,
                            const PhysicalScalarField& pressure, double dpdx, double dpdz) {
    accumulate(reduce(velocity, pressure, true, false), dpdx, dpdz);
}

void Statistics::accumulate(const std::vector<PlaneMoments>& planes, double dpdx, double dpdz) {
    if (planes.size() != static_cast<std::size_t>(m_.NYP + 1))
        throw std::runtime_error("plane moments size mismatch");
    dpdx_mean_ = update(dpdx_mean_, dpdx, count_);
    dpdz_mean_ = update(dpdz_mean_, dpdz, count_);
    const double scale = 1.0 / static_cast<double>(m_.NXTNZ);
    cpu_device().for_each(0, m_.NYP + 1, [&](int iy) {
        const auto& sum = planes[iy].sum;
        pressure_mean_[iy] = update(pressure_mean_[iy], sum[3] * scale, count_);
        pressure_second_[iy] = update(pressure_second_[iy], sum[13] * scale, count_);
        for (int c = 0; c < 3; ++c) {
            const auto i3 = static_cast<std::size_t>(c) * (m_.NYP + 1) + iy;
            velocity_mean_[i3] = update(velocity_mean_[i3], sum[c] * scale, count_);
            velocity_second_[i3] = update(velocity_second_[i3], sum[4+c] * scale, count_);
            const auto cross = static_cast<std::size_t>(c + 3) * (m_.NYP + 1) + iy;
            velocity_second_[cross] = update(velocity_second_[cross], sum[7+c] * scale, count_);
            pressure_velocity_[i3] = update(pressure_velocity_[i3], sum[10+c] * scale, count_);
            velocity_third_[i3] = update(velocity_third_[i3], sum[14+c] * scale, count_);
            velocity_fourth_[i3] = update(velocity_fourth_[i3], sum[17+c] * scale, count_);
        }
    });
    ++count_;
}

void Statistics::write(const std::filesystem::path &run_dir) const {
    for (int directory = 0; directory <= 4; ++directory)
        std::filesystem::create_directories(run_dir / "Result" / std::to_string(directory));

    if (m_.owns_zero_x_mode()) {
        const auto p0 = run_dir / "Result/0/DPDX_DPDZ.dat";
        atomic_text(p0, [&](std::ostream &out) {
            out << count_ << '\n'
                << std::scientific << std::setprecision(17) << dpdx_mean_ << ' ' << dpdz_mean_
                << '\n';
        });
    }

    const auto p1 = run_dir / "Result/1" / rank_filename("V_Pre_mean_Field_", m_.my_id, ".dat");
    const auto p2 = run_dir / "Result/2" / rank_filename("V_Pre_moment2_Field_", m_.my_id, ".dat");
    const auto p3 = run_dir / "Result/3" / rank_filename("V_Pre_moment3_Field_", m_.my_id, ".dat");
    const auto p4 = run_dir / "Result/4" / rank_filename("V_Pre_moment4_Field_", m_.my_id, ".dat");
    atomic_text(p1, [&](std::ostream &out) {
        out << std::fixed << std::setprecision(8);
        for (int iy = 0; iy <= m_.NYP; ++iy) {

            out << std::setw(19) << m_.YC[iy];
            for (int c = 0; c < 3; ++c)
                out << std::setw(19) << velocity_mean_[c * (m_.NYP + 1) + iy];
            out << std::setw(19) << pressure_mean_[iy] << '\n';
        }
    });
    atomic_text(p2, [&](std::ostream &out) {
        out << std::fixed << std::setprecision(8);
        for (int iy = 0; iy <= m_.NYP; ++iy) {

            out << std::setw(19) << m_.YC[iy];
            for (int c = 0; c < 6; ++c)
                out << std::setw(19) << velocity_second_[c * (m_.NYP + 1) + iy];
            out << std::setw(19) << pressure_second_[iy];
            for (int c = 0; c < 3; ++c)
                out << std::setw(19) << pressure_velocity_[c * (m_.NYP + 1) + iy];
            out << '\n';
        }
    });
    atomic_text(p3, [&](std::ostream &out) {
        out << std::fixed << std::setprecision(8);
        for (int iy = 0; iy <= m_.NYP; ++iy) {

            out << std::setw(19) << m_.YC[iy];
            for (int c = 0; c < 3; ++c)
                out << std::setw(19) << velocity_third_[c * (m_.NYP + 1) + iy];
            out << '\n';
        }
    });
    atomic_text(p4, [&](std::ostream &out) {
        out << std::fixed << std::setprecision(8);
        for (int iy = 0; iy <= m_.NYP; ++iy) {

            out << std::setw(19) << m_.YC[iy];
            for (int c = 0; c < 3; ++c)
                out << std::setw(19) << velocity_fourth_[c * (m_.NYP + 1) + iy];
            out << '\n';
        }
    });
}
