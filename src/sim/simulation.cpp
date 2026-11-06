#include "sim/simulation.hpp"
#include "core/cpu_device.hpp"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <complex>
#include <charconv>
#include <string>
#include <system_error>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "core/mesh.hpp"
#include "core/mpi_context.hpp"
#include "core/params.hpp"
#include "sim/initial_bc.hpp"

namespace {
std::string diagnostic_number(double value) {
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    if (result.ec != std::errc{})
        throw std::runtime_error("failed to format diagnostic value");
    return std::string(buffer, result.ptr);
}

initial_bc::DisturbanceEngine::result_type disturbance_seed(const Mesh &mesh) {
    return static_cast<initial_bc::DisturbanceEngine::result_type>(
               std::chrono::high_resolution_clock::now().time_since_epoch().count()) +
           static_cast<initial_bc::DisturbanceEngine::result_type>(1000 * mesh.my_id);
}
} // namespace

StepCoefficients step_coefficients(int iteration, double dt, double reynolds) {
    if (iteration < 1 || dt <= 0.0)
        throw std::invalid_argument("invalid time-step coefficients");
    const int effective = iteration < 3 ? iteration : 2;
    const double a0 = static_cast<double>(effective) / dt;
    const double a1 = static_cast<double>(1 - effective) / (2.0 * dt);
    const double b0 = static_cast<double>(effective);
    const double b1 = static_cast<double>(1 - effective);
    return {a0, a1, b0, b1, (a0 + a1) * reynolds};
}

Simulation::Simulation(const Params &params, const Mesh &mesh, Transforms &transforms,
                       const std::filesystem::path &run_dir, const MpiContext *mpi)
    : p_(params), m_(mesh), mpi_(mpi), transforms_(transforms), run_dir_(run_dir),
      pressure_(make_scalar_field(mesh)), velocity_(make_vector_field(mesh)),
      pressure_bc_(make_scalar_boundary(mesh)), velocity_bc_(make_vector_boundary(mesh)),
      nonlinear_(make_vector_field(mesh)), velocity_prev_(make_vector_field(mesh)),
      nonlinear_prev_(make_vector_field(mesh)), pressure_bc_prev_(make_scalar_boundary(mesh)),
      disturbance_engine_(disturbance_seed(mesh)), solver_(params, mesh, transforms),
      timestepper_(params, mesh, transforms, solver_),
      field_io_(mesh, transforms, run_dir, mpi, params.output_format), statistics_(mesh),
      diagnostics_(mesh, transforms, mpi) {
    disturbances_remaining_ = p_.ctrl.NCOUNT_disturb;
    pressure_.fill({});
    velocity_.fill({});
    pressure_bc_.fill({});
    velocity_bc_.fill({});
    nonlinear_.fill({});
    velocity_prev_.fill({});
    nonlinear_prev_.fill({});
    pressure_bc_prev_.fill({});
}

void Simulation::initialize() {
    if (p_.ctrl.INDEX_CONT == 0) {
        nt_ = 0;
        time_ = 0.0;
        initial_bc::initialize_velocity(velocity_, p_, m_, transforms_, disturbance_engine_);
    } else if (p_.ctrl.INDEX_CONT == 1) {
        field_io_.read_restart(nt_, time_, velocity_);
        if (p_.ctrl.NCOUNT_disturb > 0) {
            const auto completed = std::min<std::int64_t>(p_.ctrl.NCOUNT_disturb,
                                                         nt_ / p_.ctrl.NDISTURBFRQ);
            disturbances_remaining_ = p_.ctrl.NCOUNT_disturb - static_cast<int>(completed);
            // Restart fields are written before the injection on the same step.
            if (nt_ > 0 && nt_ % p_.ctrl.NDISTURBFRQ == 0 &&
                nt_ / p_.ctrl.NDISTURBFRQ <= p_.ctrl.NCOUNT_disturb) {
                ++disturbances_remaining_;
                apply_disturbance(nt_);
            }
        }
    } else {
        throw std::runtime_error("INDEX_CONT must be 0 or 1");
    }
    initial_bc::velocity_boundary(velocity_bc_, time_, p_, m_);
}

void Simulation::apply_disturbance(std::int64_t step) {
    initial_bc::disturb_velocity(velocity_, p_.ctrl.EPS_disturb, p_, m_, transforms_,
                                 disturbance_engine_);
    --disturbances_remaining_;
    if (m_.owns_zero_x_mode()) {
        std::cout << "disturbance at step " << step << " ("
                  << disturbances_remaining_ << " remaining)\n";
        std::cout.flush();
    }
}

double Simulation::flow_rate(const VectorField &field, int component) const {
    double flow = 0.0;
    for (int iy = 0; iy <= m_.NYS; iy += 2)
        flow -= field.at_z(0, iy, 0, component).real() / static_cast<double>(iy * iy - 1);
    return flow;
}

void Simulation::correct_flow() {
    if (m_.owns_zero_x_mode()) {
        const double qx = flow_rate(velocity_, 0);
        const double qz = flow_rate(velocity_, 2);
        if (p_.parm.INDEX_GRAD_PRE == 1) {
            dpdx_ = p_.GRAD_PRE_X;
            dpdz_ = p_.GRAD_PRE_Z;
        } else {
            const double qx_response = solver_.qx_dpdx_unit();
            const double qz_response = solver_.qz_dpdz_unit();
            if (!std::isfinite(qx_response) || qx_response == 0.0 || !std::isfinite(qz_response) ||
                qz_response == 0.0)
                throw std::runtime_error(
                    "unit pressure-gradient flow response is zero or non-finite");
            dpdx_ = (p_.QXmean - qx) / qx_response;
            dpdz_ = (p_.QZmean - qz) / qz_response;
        }
        const auto &u_response = solver_.u_dpdx_unit();
        const auto &w_response = solver_.w_dpdz_unit();
        for (int iy = 0; iy <= m_.NYS; ++iy) {
            velocity_.at_z(0, iy, 0, 0) += dpdx_ * u_response[iy];
            velocity_.at_z(0, iy, 0, 2) += dpdz_ * w_response[iy];
        }
    }
    double gradients[]{dpdx_, dpdz_};
    if (mpi_ && mpi_->size() > 1)
        mpi_->broadcast_doubles(gradients, 2, m_.zero_mode_owner(), "flow correction");
    dpdx_ = gradients[0];
    dpdz_ = gradients[1];
}

void Simulation::enqueue_serial_step(const StepCoefficients &coefficients, double step_time) {
    auto &device = cpu_device();
    timestepper_.enqueue_nonlinear_term(nonlinear_, velocity_, pressure_bc_);
    device.enqueue_for_each(-m_.NZS, m_.NZS + 1, [this, coefficients](int iz) {
            for (int ix = 0; ix <= m_.NXSP; ++ix)
                for (int wall = 0; wall < 2; ++wall)
                    pressure_bc_prev_.at_z(iz, wall, ix, 0) =
                        coefficients.b0 * pressure_bc_.at_z(iz, wall, ix, 0) +
                        coefficients.b1 * pressure_bc_prev_.at_z(iz, wall, ix, 0);
            for (int iy = 0; iy <= m_.NYS; ++iy)
                for (int ix = 0; ix <= m_.NXSP; ++ix)
                    for (int component = 0; component < 3; ++component)
                        nonlinear_prev_.at_z(iz, iy, ix, component) =
                            coefficients.a0 * velocity_.at_z(iz, iy, ix, component) +
                            coefficients.a1 * velocity_prev_.at_z(iz, iy, ix, component) +
                            coefficients.b0 * nonlinear_.at_z(iz, iy, ix, component) +
                            coefficients.b1 * nonlinear_prev_.at_z(iz, iy, ix, component);
    });
    device.enqueue([this, step_time] {
        std::swap(velocity_prev_, velocity_);
        if (p_.parm.INDEX_UNSTEADY_BC == 1)
            initial_bc::velocity_boundary(velocity_bc_, step_time, p_, m_);
    });
    timestepper_.enqueue_step(pressure_, pressure_bc_prev_, velocity_, nonlinear_prev_, velocity_bc_);
    device.enqueue([this] {
        std::swap(nonlinear_prev_, nonlinear_);
        std::swap(pressure_bc_prev_, pressure_bc_);
        correct_flow(); // rank=1 has no MPI call here.
    });
}

void Simulation::run() {
    for (int it = 1; it <= p_.ctrl.NSTEPS; ++it) {
        const auto coefficients = step_coefficients(it, p_.ctrl.DT, p_.parm.Re);
        if (it < 3)
            solver_.setup(coefficients.velocity_g);

        const std::int64_t ntpit = nt_ + it;
        time_ += p_.ctrl.DT;

        if (m_.INDEX_Serial == 1) {
            try {
                enqueue_serial_step(coefficients, time_);
                cpu_device().sync();
            } catch (...) {
                // Keep member buffers alive until already submitted work has drained.
                try { cpu_device().sync(); } catch (...) {}
                throw;
            }
        } else {
            timestepper_.nonlinear_term(nonlinear_, velocity_, pressure_bc_);
            cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
                for (int ix = 0; ix <= m_.NXSP; ++ix)
                    for (int wall = 0; wall < 2; ++wall)
                        pressure_bc_prev_.at_z(iz, wall, ix, 0) =
                            coefficients.b0 * pressure_bc_.at_z(iz, wall, ix, 0) +
                            coefficients.b1 * pressure_bc_prev_.at_z(iz, wall, ix, 0);
                for (int iy = 0; iy <= m_.NYS; ++iy)
                    for (int ix = 0; ix <= m_.NXSP; ++ix)
                        for (int component = 0; component < 3; ++component)
                            nonlinear_prev_.at_z(iz, iy, ix, component) =
                                coefficients.a0 * velocity_.at_z(iz, iy, ix, component) +
                                coefficients.a1 * velocity_prev_.at_z(iz, iy, ix, component) +
                                coefficients.b0 * nonlinear_.at_z(iz, iy, ix, component) +
                                coefficients.b1 * nonlinear_prev_.at_z(iz, iy, ix, component);
            });
            // The RHS has consumed the old history. The step overwrites every velocity
            // coefficient, so recycle its old buffer while retaining the current state.
            std::swap(velocity_prev_, velocity_);
            if (p_.parm.INDEX_UNSTEADY_BC == 1)
                initial_bc::velocity_boundary(velocity_bc_, time_, p_, m_);
            timestepper_.step(pressure_, pressure_bc_prev_, velocity_, nonlinear_prev_, velocity_bc_);
            // The next nonlinear evaluation overwrites both scratch fields in full.
            std::swap(nonlinear_prev_, nonlinear_);
            std::swap(pressure_bc_prev_, pressure_bc_);
            correct_flow();
        }

        const int physical_interval = p_.physical_interval < 0 ? p_.ctrl.NPRNFRQ : p_.physical_interval;
        const int statistics_interval = p_.statistics_interval < 0 ? p_.ctrl.NPRNFRQ : p_.statistics_interval;
        const bool diagnostic_step = ntpit % p_.ctrl.NPRNFRQ == 0;
        const bool physical_step = physical_interval > 0 && ntpit % physical_interval == 0;
        const bool restart_step = ntpit % p_.ctrl.NDUMFRQ == 0;
        const bool check_step = diagnostic_step || restart_step || it == p_.ctrl.NSTEPS;
        const bool sample = ntpit > p_.ctrl.Nstat;
        if (check_step || physical_step || sample) {
            bool finite = true;
            const auto& snapshot = field_io_.to_physical(velocity_, pressure_, check_step ? &finite : nullptr);
            const auto& planes = statistics_.reduce(snapshot.velocity, snapshot.pressure, sample, check_step);
            if (check_step) {
                const auto values = diagnostics_.evaluate(planes, velocity_, finite);
                if (diagnostic_step) diagnostics_.append(run_dir_, ntpit, time_, values);
                if (m_.owns_zero_x_mode()) {
                    std::cout << "step = " << ntpit << " TIME = " << diagnostic_number(time_)
                              << " EKT = " << diagnostic_number(values.ekt)
                              << " TKE = " << diagnostic_number(values.tke)
                              << " div_L2 = " << diagnostic_number(values.divergence_l2)
                              << " finite = " << (values.finite ? "true" : "false") << '\n';
                    std::cout.flush();
                }
                if (!values.finite)
                    throw std::runtime_error("non-finite diagnostic or field at step " + std::to_string(ntpit));
            }
            if (physical_step) field_io_.write_physical(ntpit, time_, snapshot);
            if (sample) {
                statistics_.accumulate(planes, dpdx_, dpdz_);
                if (ntpit % statistics_interval == 0 || it == p_.ctrl.NSTEPS)
                    statistics_.write(run_dir_);
            }
        }
        // Diagnostics run before overwriting a checkpoint on the same step.
        if (restart_step) field_io_.write_restart(ntpit, time_, velocity_);
        if (disturbances_remaining_ > 0 && ntpit % p_.ctrl.NDISTURBFRQ == 0)
            apply_disturbance(ntpit);
    }
    nt_ += p_.ctrl.NSTEPS;
}
