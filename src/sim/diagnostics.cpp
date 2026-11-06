#include "sim/diagnostics.hpp"
#include "core/rank_filename.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "core/mpi_context.hpp"
#include "fft/transforms.hpp"
#include "solver/vector_ops.hpp"

Diagnostics::Diagnostics(const Mesh& mesh, Transforms& transforms, const MpiContext* mpi)
    : m_(mesh), t_(transforms), mpi_(mpi), weights_(mesh.NYP + 1),
      divergence_(make_scalar_field(mesh)),
      physical_divergence_(static_cast<std::size_t>(mesh.NXTNZ) * (mesh.NYP + 1)) {
    // Clenshaw-Curtis quadrature on [-1,1]; an arithmetic y average would
    // over-weight the clustered wall points. Store weights normalized by height 2.
    const int n = mesh.NY;
    const double pi = std::acos(-1.0);
    for (int local = 0; local <= mesh.NYP; ++local) {
        const int j = mesh.y_begin(mesh.INDEX_MPI0) + local;
        const double theta = pi * j / n;
        double weight;
        if (j == 0 || j == n) weight = 1.0 / (n*n - (n % 2 == 0 ? 1.0 : 0.0));
        else {
            double value = 1.0;
            for (int k = 1; k <= (n-1)/2; ++k)
                value -= 2.0 * std::cos(2*k*theta) / (4.0*k*k-1.0);
            if (n % 2 == 0) value -= std::cos(n*theta) / (n*n-1.0);
            weight = 2.0 * value / n;
        }
        weights_[local] = weight / 2.0;
    }
}

DiagnosticValues Diagnostics::evaluate(const std::vector<PlaneMoments>& planes,
                                      const VectorField& velocity, bool spectral_finite) {
    if (planes.size() != weights_.size()) throw std::runtime_error("diagnostic plane size mismatch");
    vector_ops::divergence(divergence_, velocity, m_, t_);
    t_.channel_xyzfft('B', physical_divergence_.data(), divergence_.data());
    DiagnosticValues result;
    result.finite = spectral_finite;
    const double inv_plane = 1.0 / m_.NXTNZ;
    double div_square = 0.0;
    for (int iy = 0; iy <= m_.NYP; ++iy) {
        const auto& plane = planes[iy];
        result.finite &= plane.finite;
        result.ekt += plane.sum[5];
        double variance = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double mean = plane.sum[c] * inv_plane;
            // Roundoff in subtracting raw moments can produce a tiny negative variance.
            const double raw_variance = plane.sum[4+c] * inv_plane - mean*mean;
            result.finite &= std::isfinite(raw_variance);
            variance += std::max(0.0, raw_variance);
        }
        result.tke += 0.5 * weights_[iy] * variance;
        for (int iz = 0; iz < m_.NZ; ++iz)
            for (int ix = 0; ix < m_.NX; ++ix) {
                const auto i = physical_index(iz, iy, ix, m_.NYP+1, m_.NX);
                const double div = physical_divergence_[i];
                result.finite &= std::isfinite(div);
                div_square += weights_[iy] * inv_plane * div*div;
            }
    }
    if (mpi_ && mpi_->size() > 1) {
        const int root = m_.zero_mode_owner();
        result.ekt = mpi_->reduce_sum(result.ekt, root, "diagnostic EKT");
        result.tke = mpi_->reduce_sum(result.tke, root, "diagnostic TKE");
        div_square = mpi_->reduce_sum(div_square, root, "diagnostic divergence");
    }
    result.divergence_l2 = std::sqrt(div_square);
    result.finite &= std::isfinite(result.ekt) && std::isfinite(result.tke) &&
                     std::isfinite(result.divergence_l2);
    if (mpi_ && mpi_->size() > 1)
        result.finite = mpi_->allreduce_min(result.finite, "diagnostic finite check");
    return result;
}

void Diagnostics::append(const std::filesystem::path& directory, std::int64_t step,
                         double time, const DiagnosticValues& values) const {
    if (!m_.owns_zero_x_mode()) return;
    const auto path = directory / "diagnostics.csv";
    const bool header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
    std::ofstream csv(path, std::ios::app);
    if (header) csv << "step,time,ekt_sum_v2,tke_plane_fluctuation_volume_mean,divergence_l2,finite\n";
    csv << std::scientific << std::setprecision(17) << step << ',' << time << ','
        << values.ekt << ',' << values.tke << ',' << values.divergence_l2 << ','
        << static_cast<int>(values.finite) << '\n';
    csv.flush();
    std::ofstream ekt(directory / rank_filename("EKT", m_.my_id, ".dat"), std::ios::app);
    ekt << std::scientific << std::setprecision(17) << time << ' ' << values.ekt << '\n';
    ekt.flush();
    if (!csv || !ekt) throw std::runtime_error("cannot append diagnostics");
}
