#include <gtest/gtest.h>

#include <algorithm>

#include "fft/transforms.hpp"
#include "solver/vector_ops.hpp"
#include "stage3_fixture.hpp"

TEST(NonlinearOps, CrossSelfAndAntisymmetry) {
    const Params p = stage3_params();
    const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    VectorField a = make_vector_field(m), b = make_vector_field(m);
    VectorField aa = make_vector_field(m), ab = make_vector_field(m), ba = make_vector_field(m);
    for (int c = 0; c < 3; ++c) {
        a.at_z(0, 0, 0, c) = c + 1.0;
        b.at_z(0, 0, 0, c) = 4.0 - c;
    }
    vector_ops::cross_product(aa, a, a, m, t);
    vector_ops::cross_product(ab, a, b, m, t);
    vector_ops::cross_product(ba, b, a, m, t);
    // Independent constant-vector result, including component ordering and sign.
    const double expected[3] = {-5.0, 10.0, -5.0};
    for (int c = 0; c < 3; ++c)
        EXPECT_NEAR(std::abs(ab.at_z(0, 0, 0, c) - expected[c]), 0.0, 1e-12);
    double self = 0.0, anti = 0.0;
    for (int i = 0; i < aa.size(); ++i) {
        self = std::max(self, std::abs(aa.data()[i]));
        anti = std::max(anti, std::abs(ab.data()[i] + ba.data()[i]));
    }
    EXPECT_LT(self, 1e-13);
    EXPECT_LT(anti, 1e-13);
}
