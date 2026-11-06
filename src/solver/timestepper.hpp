#pragma once

#include "solver/eqsolve_pv.hpp"
#include "solver/vector_ops.hpp"

class TimeStepper {
  public:
    TimeStepper(const Params &params, const Mesh &mesh, Transforms &transforms, EqSolvePV &solver);

    void nonlinear_term(VectorField &nl_velocity, const VectorField &velocity,
                        ScalarBoundary &pressure_bc);
    void step(ScalarField &pressure, const ScalarBoundary &pressure_bc, VectorField &velocity,
              const VectorField &rhs_velocity, const VectorBoundary &velocity_bc);
    // Serial layout only; submitted stages are ordered by the pool queue.
    void enqueue_nonlinear_term(VectorField &nl_velocity, const VectorField &velocity,
                                ScalarBoundary &pressure_bc);
    void enqueue_step(ScalarField &pressure, const ScalarBoundary &pressure_bc, VectorField &velocity,
                       const VectorField &rhs_velocity, const VectorBoundary &velocity_bc);

  private:
    void enqueue_pressure_boundary(const VectorField &nl_velocity, ScalarBoundary &pressure_bc);
    void correct_zero_mode(VectorField &velocity, const VectorBoundary &velocity_bc);
    const Params &p_;
    const Mesh &m_;
    Transforms &t_;
    EqSolvePV &solver_;
    VectorField curl_workspace_;
    vector_ops::ProductWorkspace product_workspace_;
};
