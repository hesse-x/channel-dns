#include <gtest/gtest.h>

#include <cmath>

#include "fft/transforms.hpp"
#include "solver/eqsolve_pv.hpp"
#include "solver/vector_ops.hpp"
#include "stage3_fixture.hpp"

struct EqSolvePVTestAccess {
    using Complex = std::complex<double>;
    static void solve_tau2(EqSolvePV &solver, ScalarField &out, const ScalarField &rhs,
                           const ScalarBoundary &bc) {
        solver.solve_tau2(out, rhs, solver.g_v_, {1.0, 1.25}, {0.1, 0.2}, bc, solver.wn_velocity_,
                          solver.ad_velocity_);
    }
    static void staged(EqSolvePV &solver, ScalarField &pressure, const ScalarField &rhs_pressure,
                       const ScalarBoundary &pressure_bc, VectorField &velocity,
                       const VectorField &rhs_velocity, const VectorBoundary &velocity_bc,
                       bool correct) {

        solver.solve_coupled_base(pressure, rhs_pressure, pressure_bc, velocity, rhs_velocity,
                                  velocity_bc);
        ScalarBoundary div = make_scalar_boundary(solver.m_);
        vector_ops::divergence_boundary(div, velocity, solver.m_, solver.t_);
        cpu_device().for_each(-solver.m_.NZS, solver.m_.NZS + 1, [&](int iz) {
            for (int ix = 0; ix <= solver.m_.NXSP; ++ix) {
                const Complex c0 =
                    div.at_z(iz, 1, ix, 0) * solver.influence_divergence_[1].at_z(iz, 0, ix, 0) -
                    div.at_z(iz, 0, ix, 0) * solver.influence_divergence_[1].at_z(iz, 1, ix, 0);
                const Complex c1 =
                    div.at_z(iz, 0, ix, 0) * solver.influence_divergence_[0].at_z(iz, 1, ix, 0) -
                    div.at_z(iz, 1, ix, 0) * solver.influence_divergence_[0].at_z(iz, 0, ix, 0);
                for (int iy = 0; iy <= solver.m_.NYS; ++iy) {
                    pressure.at_z(iz, iy, ix, 0) +=
                        c0 * solver.pressure_influence_[0].at_z(iz, iy, ix, 0) +
                        c1 * solver.pressure_influence_[1].at_z(iz, iy, ix, 0);
                    for (int c = 0; c < 3; ++c)
                        velocity.at_z(iz, iy, ix, c) +=
                            c0 * solver.velocity_influence_[0].at_z(iz, iy, ix, c) +
                            c1 * solver.velocity_influence_[1].at_z(iz, iy, ix, c);
                }
            }
        });
        if (correct)
            solver.apply_sigma_correction(pressure, velocity, rhs_velocity);
    }
    static void fused(EqSolvePV &solver, ScalarField &pressure, const ScalarField &rhs_pressure,
                      const ScalarBoundary &pressure_bc, VectorField &velocity,
                      const VectorField &rhs_velocity, const VectorBoundary &velocity_bc,
                      bool correct, bool fuse_rhs = false) {
        solver.solve_modes(pressure, fuse_rhs ? nullptr : &rhs_pressure, pressure_bc, velocity,
                           rhs_velocity, velocity_bc, correct);
    }
};

TEST(EqSolve, SetupBuildsFiniteCoefficientsAndResponses) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t);
    solver.setup(150.0);
    for (int iz = -m.NZS; iz <= m.NZS; ++iz)
        for (int ix = 0; ix <= m.NXSP; ++ix) {
            EXPECT_DOUBLE_EQ(
                solver.velocity_wn().at_z(iz, 0, ix, 0),
                -(t.WAVX[ix] * t.WAVX[ix] + t.WAVZ[iz + t.zoff()] * t.WAVZ[iz + t.zoff()] + 150.0));
            for (int iy = 2; iy <= m.NYS; ++iy)
                EXPECT_TRUE(std::isfinite(solver.velocity_ad().at_z(iz, iy, ix, 0)));
        }
    EXPECT_TRUE(std::isfinite(solver.qx_dpdx_unit()));
    EXPECT_TRUE(std::isfinite(solver.qz_dpdz_unit()));
    EXPECT_NE(solver.qx_dpdx_unit(), 0.0);
    EXPECT_NE(solver.qz_dpdz_unit(), 0.0);
}

TEST(EqSolve, ZeroInputProducesFiniteConstrainedSolution) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t);
    solver.setup(150.0);
    ScalarField pressure = make_scalar_field(m), rhs_p = make_scalar_field(m);
    ScalarBoundary bc_p = make_scalar_boundary(m);
    VectorField velocity = make_vector_field(m), rhs_v = make_vector_field(m);
    VectorBoundary bc_v = make_vector_boundary(m);
    solver.solve_coupled(pressure, rhs_p, bc_p, velocity, rhs_v, bc_v);
    for (int i = 0; i < pressure.size(); ++i)
        EXPECT_TRUE(std::isfinite(std::abs(pressure.data()[i])));
    ScalarBoundary div = make_scalar_boundary(m);
    vector_ops::divergence_boundary(div, velocity, m, t);
    for (int i = 0; i < div.size(); ++i)
        EXPECT_LT(std::abs(div.data()[i]), 1e-10);
}

TEST(EqSolve, CoupledTauBoundarySolveStaysFiniteForComplexData) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t);
    solver.setup(150.0);
    ScalarField out = make_scalar_field(m), rhs = make_scalar_field(m);
    ScalarBoundary bc = make_scalar_boundary(m);
    rhs.at_z(1, 2, 1, 0) = {0.5, -0.75};
    bc.at_z(1, 0, 1, 0) = {0.2, 0.4};
    bc.at_z(1, 1, 1, 0) = {-0.3, 0.1};
    EqSolvePVTestAccess::solve_tau2(solver, out, rhs, bc);
    for (int i = 0; i < out.size(); ++i)
        EXPECT_TRUE(std::isfinite(out.data()[i].real()) && std::isfinite(out.data()[i].imag()));
}

TEST(EqSolve, SetupCanBeReusedAndUnitResponsesSatisfyNoSlip) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t);
    solver.setup(150.);
    const auto first = solver.u_dpdx_unit();
    solver.setup(300.);
    EXPECT_EQ(solver.velocity_g(), 300.);
    EXPECT_NE(solver.u_dpdx_unit(), first);
    for (const auto *response : {&solver.u_dpdx_unit(), &solver.w_dpdz_unit()}) {
        ASSERT_EQ(response->size(), static_cast<std::size_t>(m.NYS + 1));
        double upper = 0, lower = 0;
        for (int y = 0; y <= m.NYS; ++y) {
            upper += (*response)[y];
            lower += (y % 2 ? -1 : 1) * (*response)[y];
        }
        EXPECT_NEAR(upper, 0., 1e-12);
        EXPECT_NEAR(lower, 0., 1e-12);
    }
    solver.setup(150.);
    ASSERT_EQ(first.size(), solver.u_dpdx_unit().size());
    for (std::size_t i = 0; i < first.size(); ++i)
        EXPECT_NEAR(solver.u_dpdx_unit()[i], first[i], 1e-15);
}

TEST(EqSolve, ModeKernelMatchesStagedSolverForComplexColumnsAndBoundaryConditions) {
    // Cover complete tiles and different tail sizes, including a zero x mode.
    for (int nxs : {7, 8, 10, 16}) for (bool neumann_velocity : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "NXS=" << nxs
                                        << " neumann_velocity=" << neumann_velocity);
        auto p = stage3_params();
        p.mesh.NXS = nxs;
        if (neumann_velocity) {
            p.parm.a_V[1] = 0.0;
            p.parm.b_V[1] = 1.0;
        }
        Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
        Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
        EqSolvePV solver(p, m, t);
        solver.setup(150.0);
        auto rhs_p = make_scalar_field(m), rhs_v = make_vector_field(m);
        auto bc_p = make_scalar_boundary(m), bc_v = make_vector_boundary(m);
        for (int i = 0; i < rhs_p.size(); ++i)
            rhs_p.data()[i] = {std::sin(i * .13) * .001, std::cos(i * .19) * .002};
        for (int i = 0; i < rhs_v.size(); ++i)
            rhs_v.data()[i] = {std::sin(i * .17) * .001, std::cos(i * .23) * .002};
        for (int i = 0; i < bc_p.size(); ++i)
            bc_p.data()[i] = {std::sin(i * .11) * .001, std::cos(i * .29) * .002};
        for (int i = 0; i < bc_v.size(); ++i)
            bc_v.data()[i] = {std::sin(i * .31) * .001, std::cos(i * .37) * .002};
        const auto old_rhs = rhs_p, old_velocity_rhs = rhs_v;
        for (bool correct : {false, true}) for (bool fuse_rhs : {false, true}) {
            auto reference_rhs = rhs_p;
            if (fuse_rhs) vector_ops::divergence(reference_rhs, rhs_v, m, t);
            auto expected_p = make_scalar_field(m), actual_p = make_scalar_field(m);
            auto expected_v = make_vector_field(m), actual_v = make_vector_field(m);
            EqSolvePVTestAccess::staged(solver, expected_p, reference_rhs, bc_p, expected_v, rhs_v, bc_v,
                                        correct);
            EqSolvePVTestAccess::fused(solver, actual_p, rhs_p, bc_p, actual_v, rhs_v, bc_v,
                                       correct, fuse_rhs);
            for (int i = 0; i < actual_p.size(); ++i)
                ASSERT_LE(std::abs(actual_p.data()[i] - expected_p.data()[i]),
                          1e-10 * std::max(1., std::abs(expected_p.data()[i])));
            for (int i = 0; i < actual_v.size(); ++i)
                ASSERT_LE(std::abs(actual_v.data()[i] - expected_v.data()[i]),
                          1e-10 * std::max(1., std::abs(expected_v.data()[i])));
        }
        for (int i = 0; i < rhs_p.size(); ++i)
            EXPECT_EQ(rhs_p.data()[i], old_rhs.data()[i]);
        for (int i = 0; i < rhs_v.size(); ++i)
            EXPECT_EQ(rhs_v.data()[i], old_velocity_rhs.data()[i]);
    }
}
