#pragma once

#include <cstdint>
#include <filesystem>

#include "io/field_io.hpp"
#include "sim/initial_bc.hpp"
#include "sim/statistics.hpp"
#include "sim/diagnostics.hpp"
#include "solver/eqsolve_pv.hpp"
#include "solver/timestepper.hpp"

class MpiContext;

struct StepCoefficients {
    double a0, a1, b0, b1, velocity_g;
};

StepCoefficients step_coefficients(int iteration, double dt, double reynolds);

class Simulation {
public:
    Simulation(const Params& params, const Mesh& mesh, Transforms& transforms,
               const std::filesystem::path& run_dir, const MpiContext* mpi = nullptr);
    void initialize();
    void run();

    const VectorField& velocity() const { return velocity_; }
    const ScalarField& pressure() const { return pressure_; }
    std::int64_t nt() const { return nt_; }
    double time() const { return time_; }

private:
    double flow_rate(const VectorField& field, int component) const;
    void correct_flow();
    void enqueue_serial_step(const StepCoefficients &coefficients, double step_time);
    void apply_disturbance(std::int64_t step);

    const Params& p_;
    const Mesh& m_;
    const MpiContext* mpi_;
    Transforms& transforms_;
    std::filesystem::path run_dir_;
    ScalarField pressure_;
    VectorField velocity_;
    ScalarBoundary pressure_bc_;
    VectorBoundary velocity_bc_;
    VectorField nonlinear_;
    VectorField velocity_prev_;
    VectorField nonlinear_prev_;
    ScalarBoundary pressure_bc_prev_;
    std::int64_t nt_ = 0;
    double time_ = 0.0;
    int disturbances_remaining_ = 0;
    double dpdx_ = 0.0, dpdz_ = 0.0;
    initial_bc::DisturbanceEngine disturbance_engine_;
    EqSolvePV solver_;
    TimeStepper timestepper_;
    FieldIO field_io_;
    Statistics statistics_;
    Diagnostics diagnostics_;
};
