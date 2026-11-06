#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/mesh.hpp"
#include "fft/transforms.hpp"
#include "sim/initial_bc.hpp"
#include "sim/physical_fields.hpp"
#include "stage3_fixture.hpp"

TEST(InitialBC, BuildsDeterministicParabolicProfileAndBoundary) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    VectorField velocity = make_vector_field(m);
    initial_bc::initialize_velocity(velocity, p, m, transforms);

    ScalarField scalar = make_scalar_field(m);
    std::vector<double> physical(static_cast<std::size_t>(m.NZ) * (m.NYP + 1) * m.NX);
    for (int component = 0; component < 3; ++component) {
        for (int iz = -m.NZS; iz <= m.NZS; ++iz)
            for (int iy = 0; iy <= m.NYS; ++iy)
                for (int ix = 0; ix <= m.NXSP; ++ix)
                    scalar.at_z(iz, iy, ix, 0) = velocity.at_z(iz, iy, ix, component);
        transforms.channel_xyzfft('B', physical.data(), scalar.data());
        for (int iz = 0; iz < m.NZ; ++iz)
            for (int iy = 0; iy <= m.NYP; ++iy)
                for (int ix = 0; ix < m.NX; ++ix) {
                    const double expected = component == 0 ? 1.5 * (1.0 - m.YC[iy] * m.YC[iy]) :
                        (component == 1 ? p.ctrl.V_Microblow0 : 0.0);
                    EXPECT_NEAR(physical[physical_index(iz, iy, ix, m.NYP + 1, m.NX)], expected, 1e-12);
                }
    }

    VectorBoundary boundary = make_vector_boundary(m);
    boundary.fill({7.0, 2.0});
    initial_bc::velocity_boundary(boundary, 0.0, p, m);
    for (int iz = -m.NZS; iz <= m.NZS; ++iz)
        for (int wall = 0; wall < 2; ++wall)
            for (int ix = 0; ix <= m.NXSP; ++ix)
                for (int component = 0; component < 3; ++component) {
                    const auto expected = iz == 0 && ix == 0 && component == 1 ?
                        std::complex<double>(p.ctrl.V_Microblow0) : std::complex<double>();
                    EXPECT_EQ(boundary.at_z(iz, wall, ix, component), expected);
                }
}

TEST(InitialBC, AppliesRandomPerturbationWhenEnabled) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    VectorField baseline = make_vector_field(m), disturbed = make_vector_field(m), repeated = make_vector_field(m);
    initial_bc::initialize_velocity(baseline, p, m, transforms);
    p.ctrl.EPS_init = 0.1;
    initial_bc::DisturbanceEngine first(12345), second(12345);
    initial_bc::initialize_velocity(disturbed, p, m, transforms, first);
    initial_bc::initialize_velocity(repeated, p, m, transforms, second);
    double difference = 0.0;
    for (int i = 0; i < baseline.size(); ++i) {
        difference = std::max(difference, std::abs(baseline.data()[i] - disturbed.data()[i]));
        EXPECT_EQ(disturbed.data()[i], repeated.data()[i]);
    }
    EXPECT_GT(difference, 1e-6);
}

TEST(InitialBC, DisturbanceGeneratorHasUniformMomentsAndSupport) {
    initial_bc::DisturbanceEngine engine(987654321);
    constexpr int samples = 200000;
    double sum = 0.0, sum2 = 0.0, lag = 0.0;
    double previous = initial_bc::centered_uniform(engine);
    double minimum = previous, maximum = previous;
    sum += previous; sum2 += previous * previous;
    for (int i = 1; i < samples; ++i) {
        const double value = initial_bc::centered_uniform(engine);
        EXPECT_GE(value, -0.5);
        EXPECT_LT(value, 0.5);
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
        sum += value;
        sum2 += value * value;
        lag += previous * value;
        previous = value;
    }
    const double mean = sum / samples;
    const double variance = sum2 / samples - mean * mean;
    EXPECT_NEAR(mean, 0.0, 2e-3);
    EXPECT_NEAR(variance, 1.0 / 12.0, 8e-4);
    EXPECT_NEAR(lag / (samples - 1), 0.0, 8e-4);
    EXPECT_LT(minimum, -0.499);
    EXPECT_GT(maximum, 0.499);
}

TEST(InitialBC, ZeroDisturbancePreservesFieldAndRandomEngine) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    auto velocity = make_vector_field(m);
    velocity.fill({3.,4.});
    initial_bc::DisturbanceEngine engine(123), untouched(123);
    initial_bc::disturb_velocity(velocity, 0., p, m, t, engine);
    for (int i = 0; i < velocity.size(); ++i) EXPECT_EQ(velocity.data()[i], std::complex<double>(3.,4.));
    EXPECT_EQ(engine(), untouched());
}

TEST(InitialBC, NonOwnerBoundaryIsClearedAtEveryTime) {
    auto p = stage3_params();
    p.mesh.INDEX_Serial = 0;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 1, 2);
    auto boundary = make_vector_boundary(m);
    for (double time : {0., 1.25}) {
        boundary.fill({99.,88.});
        initial_bc::velocity_boundary(boundary, time, p, m);
        for (int i = 0; i < boundary.size(); ++i) EXPECT_EQ(boundary.data()[i], std::complex<double>{});
    }
}
