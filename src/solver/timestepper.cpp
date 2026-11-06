#include "solver/timestepper.hpp"
#include "core/cpu_device.hpp"

#include "solver/vector_ops.hpp"

TimeStepper::TimeStepper(const Params &params, const Mesh &mesh, Transforms &transforms,
                         EqSolvePV &solver)
    : p_(params), m_(mesh), t_(transforms), solver_(solver),
      curl_workspace_(make_vector_field(mesh)) {}

void TimeStepper::nonlinear_term(VectorField &nl_velocity, const VectorField &velocity,
                                 ScalarBoundary &pressure_bc) {
    vector_ops::curl(curl_workspace_, velocity, m_, t_);
    vector_ops::cross_product(nl_velocity, velocity, curl_workspace_, m_, t_, &product_workspace_);
    enqueue_pressure_boundary(nl_velocity, pressure_bc);
    cpu_device().sync();
}

void TimeStepper::enqueue_nonlinear_term(VectorField &nl_velocity, const VectorField &velocity,
                                         ScalarBoundary &pressure_bc) {
    vector_ops::enqueue_curl(curl_workspace_, velocity, m_, t_);
    vector_ops::enqueue_cross_product(nl_velocity, velocity, curl_workspace_, m_, t_, product_workspace_);
    enqueue_pressure_boundary(nl_velocity, pressure_bc);
}

void TimeStepper::enqueue_pressure_boundary(const VectorField &nl_velocity,
                                            ScalarBoundary &pressure_bc) {
    cpu_device().enqueue_for_each(-m_.NZS, m_.NZS + 1, [this, &nl_velocity, &pressure_bc](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            pressure_bc.at_z(iz, 0, ix, 0) = 0.0;
            pressure_bc.at_z(iz, 1, ix, 0) = 0.0;
        }
        for (int iy = 0; iy <= m_.NYS; ++iy)
            for (int ix = 0; ix <= m_.NXSP; ++ix) {
                const auto wall_term =
                    t_.CWAVXRe[ix] * curl_workspace_.at_z(iz, iy, ix, 2) -
                    t_.CWAVZRe[iz + t_.zoff()] * curl_workspace_.at_z(iz, iy, ix, 0) +
                    nl_velocity.at_z(iz, iy, ix, 1);
                pressure_bc.at_z(iz, 0, ix, 0) += wall_term;
                pressure_bc.at_z(iz, 1, ix, 0) += (iy % 2 == 0) ? -wall_term : wall_term;
            }
    });
}

void TimeStepper::step(ScalarField &pressure, const ScalarBoundary &pressure_bc,
                       VectorField &velocity, const VectorField &rhs_velocity,
                       const VectorBoundary &velocity_bc) {
    solver_.solve_modes(pressure, nullptr, pressure_bc, velocity, rhs_velocity,
                        velocity_bc, true);
    correct_zero_mode(velocity, velocity_bc);
}

void TimeStepper::enqueue_step(ScalarField &pressure, const ScalarBoundary &pressure_bc,
                                VectorField &velocity, const VectorField &rhs_velocity,
                                const VectorBoundary &velocity_bc) {
    solver_.solve_modes(pressure, nullptr, pressure_bc, velocity, rhs_velocity, velocity_bc, true, false);
    cpu_device().enqueue([this, &velocity, &velocity_bc] { correct_zero_mode(velocity, velocity_bc); });
}

void TimeStepper::correct_zero_mode(VectorField &velocity, const VectorBoundary &velocity_bc) {
    if (m_.owns_zero_x_mode()) {
        const auto lower = velocity_bc.at_z(0, 0, 0, 1);
        const auto upper = velocity_bc.at_z(0, 1, 0, 1);
        velocity.at_z(0, 0, 0, 1) = 0.5 * (lower + upper);
        for (int iy = 1; iy <= m_.NYS; ++iy)
            velocity.at_z(0, iy, 0, 1) = 0.0;
    }
}
