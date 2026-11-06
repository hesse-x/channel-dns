#include <gtest/gtest.h>

#include "fft/transforms.hpp"
#include "solver/timestepper.hpp"
#include "stage3_fixture.hpp"

TEST(TimeStepper, ZeroNonlinearTermIsRepeatable) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t); solver.setup(150.0);
    TimeStepper stepper(p, m, t, solver);
    VectorField velocity = make_vector_field(m), nl = make_vector_field(m);
    ScalarBoundary bc = make_scalar_boundary(m);
    nl.fill({3.0, 4.0}); bc.fill({3.0, 4.0});
    stepper.nonlinear_term(nl, velocity, bc);
    for (int i = 0; i < nl.size(); ++i) EXPECT_EQ(nl.data()[i], (std::complex<double>{}));
    for (int i = 0; i < bc.size(); ++i) EXPECT_EQ(bc.data()[i], (std::complex<double>{}));
}

TEST(TimeStepper, ZeroStepStaysFinite) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t); solver.setup(150.0);
    TimeStepper stepper(p, m, t, solver);
    ScalarField pressure = make_scalar_field(m);
    ScalarBoundary pressure_bc = make_scalar_boundary(m);
    VectorField velocity = make_vector_field(m), rhs = make_vector_field(m);
    VectorBoundary velocity_bc = make_vector_boundary(m);
    stepper.step(pressure, pressure_bc, velocity, rhs, velocity_bc);
    for (int i = 0; i < pressure.size(); ++i) EXPECT_TRUE(std::isfinite(std::abs(pressure.data()[i])));
    for (int i = 0; i < velocity.size(); ++i) EXPECT_TRUE(std::isfinite(std::abs(velocity.data()[i])));
}

TEST(TimeStepper, ProjectsMeanWallNormalVelocityToConstantMode) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p, m, t); solver.setup(150.0);
    TimeStepper stepper(p, m, t, solver);
    ScalarField pressure = make_scalar_field(m);
    ScalarBoundary pressure_bc = make_scalar_boundary(m);
    VectorField velocity = make_vector_field(m), rhs = make_vector_field(m);
    VectorBoundary velocity_bc = make_vector_boundary(m);
    velocity_bc.at_z(0, 0, 0, 1) = 1.0e-3;
    velocity_bc.at_z(0, 1, 0, 1) = 1.0e-3;

    stepper.step(pressure, pressure_bc, velocity, rhs, velocity_bc);

    EXPECT_EQ(velocity.at_z(0, 0, 0, 1), std::complex<double>(1.0e-3, 0.0));
    for (int iy = 1; iy <= m.NYS; ++iy)
        EXPECT_EQ(velocity.at_z(0, iy, 0, 1), std::complex<double>());
}

TEST(TimeStepper, RepeatedNonlinearCallsOverwriteWorkspaceAndPreserveVelocity) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    EqSolvePV solver(p,m,t); solver.setup(150.);
    TimeStepper stepper(p,m,t,solver);
    auto velocity = make_vector_field(m), nl = make_vector_field(m);
    auto bc = make_scalar_boundary(m);
    // Uniform velocity has zero curl and nonlinear term.
    for (int c = 0; c < 3; ++c) velocity.at_z(0,0,0,c) = c+1.;
    const auto original = velocity;
    for (int attempt = 0; attempt < 2; ++attempt) {
        nl.fill({99.,88.}); bc.fill({99.,88.});
        stepper.nonlinear_term(nl, velocity, bc);
        for (int i = 0; i < nl.size(); ++i) EXPECT_LT(std::abs(nl.data()[i]), 1e-12);
        for (int i = 0; i < bc.size(); ++i) EXPECT_LT(std::abs(bc.data()[i]), 1e-12);
        for (int i = 0; i < velocity.size(); ++i) EXPECT_EQ(velocity.data()[i], original.data()[i]);
    }
}
