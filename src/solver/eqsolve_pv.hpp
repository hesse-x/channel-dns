#pragma once

#include <array>
#include <vector>

#include "core/params.hpp"
#include "solver/fields.hpp"

class Transforms;
class TimeStepper;
struct EqSolvePVTestAccess;

class EqSolvePV {
  public:
    EqSolvePV(const Params &params, const Mesh &mesh, const Transforms &transforms);

    void setup(double g_v);
    void solve_coupled(ScalarField &pressure, const ScalarField &rhs_pressure,
                       const ScalarBoundary &pressure_bc, VectorField &velocity,
                       const VectorField &rhs_velocity, const VectorBoundary &velocity_bc);

    double velocity_g() const { return g_v_; }
    const RealModeField &pressure_wn() const { return wn_pressure_; }
    const RealModeField &velocity_wn() const { return wn_velocity_; }
    const RealScalarField &pressure_ad() const { return ad_pressure_; }
    const RealScalarField &velocity_ad() const { return ad_velocity_; }
    const std::vector<double> &u_dpdx_unit() const { return u_dpdx_unit_; }
    const std::vector<double> &w_dpdz_unit() const { return w_dpdz_unit_; }
    double qx_dpdx_unit() const { return qx_dpdx_unit_; }
    double qz_dpdz_unit() const { return qz_dpdz_unit_; }

  private:
    friend class TimeStepper;
    friend struct EqSolvePVTestAccess;
    using Pair = std::array<double, 2>;

    void build_wn_ad(RealModeField &wn, RealScalarField &ad, double g);
    void solve_tau(ScalarField &out, const ScalarField &rhs, double g, double a, double b,
                   const ScalarBoundary &bc, const RealModeField &wn, const RealScalarField &ad,
                   int component = 0) const;
    void solve_tau2(ScalarField &out, const ScalarField &rhs, double g, const Pair &a,
                    const Pair &b, const ScalarBoundary &bc, const RealModeField &wn,
                    const RealScalarField &ad) const;
    void solve_coupled_base(ScalarField &pressure, const ScalarField &rhs_pressure,
                            const ScalarBoundary &pressure_bc, VectorField &velocity,
                            const VectorField &rhs_velocity, const VectorBoundary &velocity_bc);
    // A null pressure RHS computes divergence locally from the cached velocity RHS.
    void solve_modes(ScalarField &pressure, const ScalarField *rhs_pressure,
                     const ScalarBoundary &pressure_bc, VectorField &velocity,
                     const VectorField &rhs_velocity, const VectorBoundary &velocity_bc,
                     bool correct_sigma, bool wait = true);
    void build_influence_coefficients();
    void build_b_problem();
    void build_pressure_gradient_responses();
    void apply_sigma_correction(ScalarField &pressure, VectorField &velocity,
                                const VectorField &rhs_velocity);

    const Params &p_;
    const Mesh &m_;
    const Transforms &t_;
    double g_v_ = 0.0;
    double pressure_a_ = 0.0, pressure_bcoef_ = 0.0;
    RealModeField wn_pressure_, wn_velocity_;
    RealScalarField ad_pressure_, ad_velocity_;
    std::array<ScalarField, 2> pressure_influence_;
    std::array<VectorField, 2> velocity_influence_;
    std::array<ScalarBoundary, 2> influence_divergence_;
    ScalarField pressure_b_;
    VectorField velocity_b_;
    ScalarBoundary sigma0_;
    VectorField gradient_workspace_;
    ScalarField forcing_workspace_;
    ScalarBoundary component_bc_workspace_;
    std::vector<double> u_dpdx_unit_, w_dpdz_unit_;
    double qx_dpdx_unit_ = 0.0, qz_dpdz_unit_ = 0.0;
};
