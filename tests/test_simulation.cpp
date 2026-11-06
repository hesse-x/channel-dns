#include <gtest/gtest.h>

#include "sim/simulation.hpp"

TEST(Simulation, StepCoefficientsFollowStartupScheme) {
    const auto first = step_coefficients(1, 0.25, 10.0);
    EXPECT_DOUBLE_EQ(first.a0, 4.0);
    EXPECT_DOUBLE_EQ(first.a1, 0.0);
    EXPECT_DOUBLE_EQ(first.b0, 1.0);
    EXPECT_DOUBLE_EQ(first.b1, 0.0);
    EXPECT_DOUBLE_EQ(first.velocity_g, 40.0);
    const auto second = step_coefficients(2, 0.25, 10.0);
    const auto third = step_coefficients(3, 0.25, 10.0);
    EXPECT_DOUBLE_EQ(second.a0, 8.0);
    EXPECT_DOUBLE_EQ(second.a1, -2.0);
    EXPECT_DOUBLE_EQ(second.b0, 2.0);
    EXPECT_DOUBLE_EQ(second.b1, -1.0);
    EXPECT_DOUBLE_EQ(third.a0, second.a0);
    EXPECT_DOUBLE_EQ(third.a1, second.a1);
    EXPECT_DOUBLE_EQ(third.velocity_g, second.velocity_g);
}

#include <cmath>
#include "core/mesh.hpp"
#include "fft/transforms.hpp"
#include "stage3_fixture.hpp"
#include "test_support.hpp"

TEST(Simulation, InvalidStepCoefficientsAreRejected) {
    EXPECT_THROW(step_coefficients(0, .1, 10), std::invalid_argument);
    EXPECT_THROW(step_coefficients(-1, .1, 10), std::invalid_argument);
    EXPECT_THROW(step_coefficients(1, 0, 10), std::invalid_argument);
    EXPECT_THROW(step_coefficients(1, -.1, 10), std::invalid_argument);
}

TEST(Simulation, InitializeRunOutputCadenceAndContinuation) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    p.ctrl.INDEX_CONT = 0;
    p.ctrl.NSTEPS = 3;
    p.ctrl.NDUMFRQ = 2;
    p.ctrl.NPRNFRQ = 2;
    p.ctrl.Nstat = 0;
    p.output_format = OutputFormat::Both;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    Simulation sim(p, m, t, directory.path);
    sim.initialize();
    EXPECT_EQ(sim.nt(), 0);
    EXPECT_EQ(sim.time(), 0);
    EXPECT_EQ(sim.velocity().size(), make_vector_field(m).size());
    EXPECT_EQ(sim.pressure().size(), make_scalar_field(m).size());
    EXPECT_FALSE(std::filesystem::exists(directory.path / "field_S_ 0.dat"));
    sim.run();
    EXPECT_EQ(sim.nt(), 3);
    EXPECT_NEAR(sim.time(), 3*p.ctrl.DT, 1e-16);
    for (auto extension : {"dat", "bin"}) {
        EXPECT_TRUE(std::filesystem::exists(directory.path / (std::string("field_S_ 0.") + extension)));
        EXPECT_TRUE(std::filesystem::exists(directory.path / (std::string("field_P_ 0.") + extension)));
    }
    auto velocity = make_vector_field(m);
    FieldIO io(m, t, directory.path);
    std::int64_t step;
    double time;
    io.read_restart(step, time, velocity);
    EXPECT_EQ(step, 2); // Last dump, not the last simulated step.
    EXPECT_NEAR(time, 2*p.ctrl.DT, 1e-16);
    std::ifstream ekt(directory.path / "EKT 0.dat");
    double et, energy, extra;
    ASSERT_TRUE(ekt >> et >> energy);
    EXPECT_EQ(et, time);
    EXPECT_FALSE(ekt >> extra);
    EXPECT_TRUE(std::filesystem::exists(directory.path / "Result/1/V_Pre_mean_Field_ 0.dat"));

    std::ifstream stats(directory.path / "Result/0/DPDX_DPDZ.dat");
    std::int64_t samples = 0;
    ASSERT_TRUE(stats >> samples);
    EXPECT_EQ(samples, 3); // Every step, including steps without field output.

    auto restart = p;
    restart.ctrl.INDEX_CONT = 1;
    restart.ctrl.NSTEPS = 1;
    Simulation resumed(restart, m, t, directory.path);
    resumed.initialize();
    EXPECT_EQ(resumed.nt(), 2);
    EXPECT_EQ(resumed.time(), time);
    for (int i = 0; i < velocity.size(); ++i) EXPECT_EQ(resumed.velocity().data()[i], velocity.data()[i]);
    resumed.run();
    EXPECT_EQ(resumed.nt(), 3);
    EXPECT_NEAR(resumed.time(), sim.time(), 1e-16);
}

TEST(Simulation, ZeroStepsDoesNotWriteFilesAndInvalidContinuationFails) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    p.ctrl.NSTEPS = 0;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    Simulation sim(p, m, t, directory.path);
    sim.initialize();
    const auto original = sim.velocity();
    sim.run();
    EXPECT_TRUE(std::filesystem::is_empty(directory.path));
    EXPECT_EQ(sim.nt(), 0);
    for (int i = 0; i < original.size(); ++i) EXPECT_EQ(sim.velocity().data()[i], original.data()[i]);
    p.ctrl.INDEX_CONT = 2;
    Simulation invalid(p, m, t, directory.path);
    EXPECT_THROW(invalid.initialize(), std::runtime_error);
    p.ctrl.INDEX_CONT = 1;
    Simulation missing(p, m, t, directory.path);
    EXPECT_THROW(missing.initialize(), std::runtime_error);
}

TEST(Simulation, ScalarDiagnosticsWithoutPhysicalDumpAndIndependentStatisticsCadence) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    p.ctrl.NSTEPS = 3;
    p.ctrl.Nstat = 0;
    p.ctrl.NPRNFRQ = 1;
    p.ctrl.NDUMFRQ = 2;
    p.physical_interval = 0;
    p.statistics_interval = 2;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    Simulation simulation(p, m, t, directory.path);
    simulation.initialize();
    simulation.run();
    EXPECT_FALSE(std::filesystem::exists(directory.path / "field_P_ 0.dat"));
    EXPECT_TRUE(std::filesystem::exists(directory.path / "field_S_ 0.dat"));
    std::ifstream csv(directory.path / "diagnostics.csv");
    std::string line;
    int rows = 0;
    while (std::getline(csv, line)) ++rows;
    EXPECT_EQ(rows, 4); // Header plus three diagnostic steps.
    std::ifstream stats(directory.path / "Result/0/DPDX_DPDZ.dat");
    std::int64_t count = 0;
    ASSERT_TRUE(stats >> count);
    EXPECT_EQ(count, 3); // Final save includes the last sample between save intervals.
}

TEST(Simulation, RestartRestoresDisturbanceScheduleUsingGlobalSteps) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    p.ctrl.NSTEPS = 2;
    p.ctrl.NDUMFRQ = 2;
    p.ctrl.NPRNFRQ = 100;
    p.ctrl.Nstat = 100;
    p.ctrl.NCOUNT_disturb = 2;
    p.ctrl.NDISTURBFRQ = 2;
    p.ctrl.EPS_disturb = .001;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    Simulation first(p,m,t,directory.path);
    first.initialize();
    first.run(); // Checkpoint 2 is written before injection 2.
    p.ctrl.INDEX_CONT = 1;
    p.ctrl.NSTEPS = 4;
    testing::internal::CaptureStdout();
    Simulation resumed(p,m,t,directory.path);
    resumed.initialize();
    resumed.run();
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_NE(log.find("disturbance at step 2 (1 remaining)"), std::string::npos);
    EXPECT_NE(log.find("disturbance at step 4 (0 remaining)"), std::string::npos);
    EXPECT_EQ(log.find("disturbance at step 6"), std::string::npos);
    EXPECT_EQ(resumed.nt(), 6);
}

TEST(Simulation, SmallGridSavedRestartMatchesContinuousRun) {
    for (const auto format : {OutputFormat::Ascii, OutputFormat::Binary}) {
        SCOPED_TRACE(format == OutputFormat::Ascii ? "ASCII" : "binary");
        test_support::TempDirectory continuous_dir, restart_dir;
        auto p = stage3_params();
        p.mesh = {8, 8, 8, 3, 4, 6, 4, 8, 8, 8, 1};
        p.ctrl.INDEX_CONT = 1;
        p.ctrl.NSTEPS = 6;
        p.ctrl.NDUMFRQ = 3;
        p.ctrl.NPRNFRQ = 100;
        p.ctrl.Nstat = 100;
        p.ctrl.NCOUNT_disturb = 0;
        p.ctrl.EPS_init = 0.01;
        p.ctrl.DT = 1e-4;
        p.physical_interval = 0;
        p.output_format = format;
        Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
        Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
        auto initial = make_vector_field(m);
        initial_bc::DisturbanceEngine engine(12345);
        initial_bc::initialize_velocity(initial, p, m, transforms, engine);
        for (const auto& dir : {continuous_dir.path, restart_dir.path}) {
            FieldIO io(m, transforms, dir, nullptr, format);
            io.write_restart(0, 0.0, initial);
        }
        Simulation continuous(p, m, transforms, continuous_dir.path);
        continuous.initialize();
        continuous.run();

        auto first_params = p;
        first_params.ctrl.NSTEPS = 3;
        auto checkpoint = make_vector_field(m);
        {
            Simulation first(first_params, m, transforms, restart_dir.path);
            first.initialize();
            first.run();
            ASSERT_EQ(first.nt(), 3);
            checkpoint = first.velocity();
        } // Destroy the first simulation: continuation must read its saved state.
        auto resumed_params = first_params;
        Simulation resumed(resumed_params, m, transforms, restart_dir.path);
        resumed.initialize();
        ASSERT_EQ(resumed.nt(), 3);
        EXPECT_NEAR(resumed.time(), 3 * p.ctrl.DT, 1e-18);
        for (int i = 0; i < checkpoint.size(); ++i)
            ASSERT_EQ(resumed.velocity().data()[i], checkpoint.data()[i]);
        resumed.run();
        EXPECT_EQ(resumed.nt(), continuous.nt());
        EXPECT_DOUBLE_EQ(resumed.time(), continuous.time());
        double error2 = 0.0, reference2 = 0.0, evolution2 = 0.0;
        for (int i = 0; i < initial.size(); ++i) {
            const auto expected = continuous.velocity().data()[i];
            const auto actual = resumed.velocity().data()[i];
            ASSERT_TRUE(std::isfinite(actual.real()) && std::isfinite(actual.imag()));
            error2 += std::norm(actual - expected);
            reference2 += std::norm(expected);
            evolution2 += std::norm(expected - initial.data()[i]);
        }
        ASSERT_GT(reference2, 0.0);
        EXPECT_GT(evolution2, 1e-16); // Exercise an evolving, perturbed field.
        // Restart files contain velocity, not multistep history, so startup is repeated.
        EXPECT_LT(std::sqrt(error2 / reference2), 1e-7);
    }
}
