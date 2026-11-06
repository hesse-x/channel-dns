#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "core/complex_ops.hpp"
#include "fft/transforms.hpp"
#include "solver/vector_ops.hpp"
#include "stage3_fixture.hpp"

TEST(VectorKernels, FourierDerivativeMatchesComplexMultiplication) {
    for (double k : {-3., -.1, 0., .1, 3.}) for (int i = 0; i < 31; ++i) {
        const std::complex<double> v{std::sin(i * .3), std::cos(i * .7)};
        EXPECT_EQ(fourier_derivative(k, v), std::complex<double>(0., k) * v);
    }
}

TEST(VectorOps, FieldShapeAndFourierDerivatives) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    ScalarField f = make_scalar_field(m), d = make_scalar_field(m);
    f.at_z(-2, 3, 2, 0) = {0.25, -0.75};
    vector_ops::derivative(0, d, f, m, t);
    EXPECT_NEAR(std::abs(d.at_z(-2, 3, 2, 0) - t.CWAVX[2] * f.at_z(-2, 3, 2, 0)), 0.0, 1e-14);
    vector_ops::derivative(2, d, f, m, t);
    EXPECT_NEAR(std::abs(d.at_z(-2, 3, 2, 0) - t.CWAVZ[-2 + t.zoff()] * f.at_z(-2, 3, 2, 0)), 0.0,
                1e-14);
    EXPECT_EQ(f.z_offset(), m.NZS);
    EXPECT_EQ(f.dim(0), m.NZS2P1);
}

TEST(VectorOps, ChebyshevDerivativeAndFlowRate) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    ScalarField f = make_scalar_field(m), d = make_scalar_field(m);
    f.at_z(0, 4, 0, 0) = 1.0;
    vector_ops::derivative(1, d, f, m, t);
    EXPECT_DOUBLE_EQ(d.at_z(0, 3, 0, 0).real(), 8.0);
    EXPECT_DOUBLE_EQ(d.at_z(0, 1, 0, 0).real(), 8.0);
    f.fill({});
    f.at_z(0, 0, 0, 0) = 2.0;
    f.at_z(0, 2, 0, 0) = 3.0;
    EXPECT_DOUBLE_EQ(vector_ops::flow_rate(f, m), 1.0);
}

TEST(VectorOps, DivCurlVanishes) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    VectorField v = make_vector_field(m), c = make_vector_field(m);
    ScalarField div = make_scalar_field(m);
    for (int iz = -2; iz <= 2; ++iz)
        for (int iy = 0; iy <= 4; ++iy)
            for (int ix = 0; ix <= m.NXSP; ++ix)
                for (int k = 0; k < 3; ++k)
                    v.at_z(iz, iy, ix, k) = {0.01 * (1 + ix + iy + k), 0.02 * iz};
    vector_ops::curl(c, v, m, t);
    vector_ops::divergence(div, c, m, t);
    ScalarField zero = make_scalar_field(m);
    EXPECT_LT(vector_ops::max_distance(div, zero, m), 1e-12);
}

TEST(VectorOps, CurlMatchesIndependentScalarDerivatives) {
    for (int nys : {4, 7, 30}) for (int nxs : {7, 8, 16}) {
        SCOPED_TRACE(::testing::Message() << "NYS=" << nys << " NXS=" << nxs);
        auto p = stage3_params();
        p.mesh.NYS = nys; p.mesh.NXS = nxs;
        const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
        Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
        auto input = make_vector_field(m), result = make_vector_field(m);
        std::array<ScalarField, 3> components;
        std::array<std::array<ScalarField, 3>, 3> derivatives;
        for (int c = 0; c < 3; ++c) {
            components[c] = make_scalar_field(m);
            for (int z = -m.NZS; z <= m.NZS; ++z)
                for (int y = 0; y <= m.NYS; ++y)
                    for (int x = 0; x <= m.NXSP; ++x) {
                        const auto value = std::complex<double>{std::sin(.3 * x + .7 * y + c),
                                                                std::cos(.2 * z + .9 * y + c)};
                        input.at_z(z, y, x, c) = components[c].at_z(z, y, x, 0) = value;
                    }
            for (int axis = 0; axis < 3; ++axis) {
                derivatives[c][axis] = make_scalar_field(m);
                vector_ops::derivative(axis, derivatives[c][axis], components[c], m, t);
            }
        }
        const auto original = input;
        vector_ops::curl(result, input, m, t);
        for (int c = 0; c < 3; ++c)
            for (int z = -m.NZS; z <= m.NZS; ++z)
                for (int y = 0; y <= m.NYS; ++y)
                    for (int x = 0; x <= m.NXSP; ++x) {
                        const int j = (c + 1) % 3, k = (c + 2) % 3;
                        const auto expected = derivatives[k][j].at_z(z, y, x, 0) -
                                              derivatives[j][k].at_z(z, y, x, 0);
                        EXPECT_NEAR(std::abs(result.at_z(z, y, x, c) - expected), 0., 1e-12);
                        EXPECT_EQ(input.at_z(z, y, x, c), original.at_z(z, y, x, c));
                    }
    }
}

TEST(VectorOps, GradientAndLaplacianMatchManufacturedFourierChebyshevMode) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    auto scalar = make_scalar_field(m), laplacian = make_scalar_field(m);
    auto gradient = make_vector_field(m);
    // f = T2(y) exp(i*kx*x+i*kz*z); T2'=4*T1 and T2''=4*T0.
    scalar.at_z(-2, 2, 1, 0) = {2., -1.};
    const auto input = scalar;
    vector_ops::gradient(gradient, scalar, m, t);
    vector_ops::laplacian(laplacian, scalar, m, t);
    for (int z = -m.NZS; z <= m.NZS; ++z)
        for (int y = 0; y <= m.NYS; ++y)
            for (int x = 0; x <= m.NXSP; ++x) {
                const bool mode = z == -2 && x == 1;
                const auto amplitude = std::complex<double>(2., -1.);
                const auto dx = mode && y == 2 ? t.CWAVX[1] * amplitude : std::complex<double>{};
                const auto dy = mode && y == 1 ? 4. * amplitude : std::complex<double>{};
                const auto dz =
                    mode && y == 2 ? t.CWAVZ[-2 + t.zoff()] * amplitude : std::complex<double>{};
                auto expected_laplacian = std::complex<double>{};
                if (mode && y == 0)
                    expected_laplacian = 4. * amplitude;
                if (mode && y == 2)
                    expected_laplacian =
                        -(t.WAVX[1] * t.WAVX[1] + t.WAVZ[-2 + t.zoff()] * t.WAVZ[-2 + t.zoff()]) *
                        amplitude;
                EXPECT_NEAR(std::abs(gradient.at_z(z, y, x, 0) - dx), 0, 1e-13);
                EXPECT_NEAR(std::abs(gradient.at_z(z, y, x, 1) - dy), 0, 1e-13);
                EXPECT_NEAR(std::abs(gradient.at_z(z, y, x, 2) - dz), 0, 1e-13);
                EXPECT_NEAR(std::abs(laplacian.at_z(z, y, x, 0) - expected_laplacian), 0, 1e-12);
            }
    for (int i = 0; i < scalar.size(); ++i)
        EXPECT_EQ(scalar.data()[i], input.data()[i]);
}

TEST(VectorOps, DotProductMatchesConstantVectorAlgebra) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    auto a = make_vector_field(m), b = make_vector_field(m);
    auto out = make_scalar_field(m);
    for (int c = 0; c < 3; ++c) {
        a.at_z(0, 0, 0, c) = c + 1.;
        b.at_z(0, 0, 0, c) = c + 4.;
    }
    const auto old_a = a, old_b = b;
    vector_ops::dot_product(out, a, b, m, t);
    EXPECT_NEAR(out.at_z(0, 0, 0, 0).real(), 32., 1e-12);
    out.at_z(0, 0, 0, 0) -= 32.;
    for (int i = 0; i < out.size(); ++i)
        EXPECT_LT(std::abs(out.data()[i]), 1e-12);
    for (int i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a.data()[i], old_a.data()[i]);
        EXPECT_EQ(b.data()[i], old_b.data()[i]);
    }
}

TEST(VectorOps, FusedKernelsMatchIndependentDerivativeComposition) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    auto v = make_vector_field(m), got_curl = make_vector_field(m);
    auto got_div = make_scalar_field(m), scalar = make_scalar_field(m);
    auto d = make_scalar_field(m), expected_div = make_scalar_field(m);
    auto expected_curl = make_vector_field(m);
    for (int i = 0; i < v.size(); ++i)
        v.data()[i] = {std::sin(i * .17), std::cos(i * .23)};
    for (int c = 0; c < 3; ++c) {
        for (int z = -m.NZS; z <= m.NZS; ++z)
            for (int y = 0; y <= m.NYS; ++y)
                for (int x = 0; x <= m.NXSP; ++x)
                    scalar.at_z(z, y, x, 0) = v.at_z(z, y, x, c);
        vector_ops::derivative(c, d, scalar, m, t);
        for (int i = 0; i < d.size(); ++i)
            expected_div.data()[i] += d.data()[i];
        for (int axis = 0; axis < 3; ++axis) {
            if (axis == c)
                continue;
            vector_ops::derivative(axis, d, scalar, m, t);
            const int component = 3 - axis - c;
            const double sign = (axis + 1) % 3 == c ? 1.0 : -1.0;
            for (int z = -m.NZS; z <= m.NZS; ++z)
                for (int y = 0; y <= m.NYS; ++y)
                    for (int x = 0; x <= m.NXSP; ++x)
                        expected_curl.at_z(z, y, x, component) += sign * d.at_z(z, y, x, 0);
        }
    }
    vector_ops::divergence(got_div, v, m, t);
    vector_ops::curl(got_curl, v, m, t);
    for (int i = 0; i < got_div.size(); ++i)
        EXPECT_NEAR(std::abs(got_div.data()[i] - expected_div.data()[i]), 0.0, 1e-11);
    for (int i = 0; i < got_curl.size(); ++i)
        EXPECT_NEAR(std::abs(got_curl.data()[i] - expected_curl.data()[i]), 0.0, 1e-11);
}
