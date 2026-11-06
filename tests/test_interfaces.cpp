#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <type_traits>

#include "fft/transforms.hpp"
#include "fft/wisdom.hpp"
#include "solver/vector_ops.hpp"
#include "stage3_fixture.hpp"
#include "test_support.hpp"

TEST(ModuleInterfaces, FieldFactoriesMatchSpectralAndBoundaryLayouts) {
    auto p = stage3_params();
    p.mesh.INDEX_Serial = 0;
    for (int rank = 0; rank < 3; ++rank) {
        const Mesh m(p.mesh, p.parm.XL, p.parm.ZL, rank, 3);
        auto check = [&](const auto& field, int ny, int components) {
            EXPECT_EQ(field.dim(0), m.NZS2P1);
            EXPECT_EQ(field.dim(1), ny);
            EXPECT_EQ(field.dim(2), m.x_count(rank));
            EXPECT_EQ(field.dim(3), components);
            EXPECT_EQ(field.z_offset(), m.NZS);
            EXPECT_EQ(field.size(), m.NZS2P1 * ny * m.x_count(rank) * components);
            EXPECT_EQ(field.raw(), field.data());
            for (int i = 0; i < field.size(); ++i) EXPECT_EQ(field.data()[i], typename std::decay_t<decltype(field.data()[i])>{});
        };
        check(make_scalar_field(m), m.NYS + 1, 1);
        check(make_vector_field(m), m.NYS + 1, 3);
        check(make_scalar_boundary(m), 2, 1);
        check(make_vector_boundary(m), 2, 3);
        check(make_real_scalar_field(m), m.NYS + 1, 1);
        check(make_real_mode_field(m), 1, 1);
    }
}

TEST(ModuleInterfaces, MeshCoordinateSlicesMatchGlobalChebyshevAndPeriodicGrid) {
    auto p = stage3_params();
    p.mesh.INDEX_Serial = 0;
    for (int rank = 0; rank < 3; ++rank) {
        const Mesh m(p.mesh, 2., 3., rank, 3);
        ASSERT_EQ(m.YC.size(), static_cast<std::size_t>(m.y_count(rank)));
        for (int i = 0; i < m.y_count(rank); ++i)
            EXPECT_NEAR(m.YC[i], std::cos(std::acos(-1.) * (m.y_begin(rank)+i) / m.NY), 1e-15);
        EXPECT_DOUBLE_EQ(m.XC.front(), 2. / m.NX);
        EXPECT_DOUBLE_EQ(m.XC.back(), 2.);
        EXPECT_DOUBLE_EQ(m.ZC.front(), 3. / m.NZ);
        EXPECT_DOUBLE_EQ(m.ZC.back(), 3.);
        for (int bad : {-1, 3}) {
            EXPECT_THROW(m.x_begin(bad), std::out_of_range);
            EXPECT_THROW(m.x_count(bad), std::out_of_range);
            EXPECT_THROW(m.y_begin(bad), std::out_of_range);
            EXPECT_THROW(m.y_count(bad), std::out_of_range);
            EXPECT_THROW(m.y2_begin(bad), std::out_of_range);
            EXPECT_THROW(m.y2_count(bad), std::out_of_range);
        }
    }
}

TEST(ModuleInterfaces, MeshRejectsInvalidRankSizeAndSerialFlag) {
    auto p = stage3_params();
    EXPECT_THROW((Mesh(p.mesh, 2., 3., -1, 1)), std::invalid_argument);
    EXPECT_THROW((Mesh(p.mesh, 2., 3., 1, 1)), std::invalid_argument);
    EXPECT_THROW((Mesh(p.mesh, 2., 3., 0, 0)), std::invalid_argument);
    p.mesh.INDEX_Serial = 2;
    EXPECT_THROW((Mesh(p.mesh, 2., 3., 0, 1)), std::invalid_argument);
}

TEST(ModuleInterfaces, DistributedTransformsRequireMatchingMpiContext) {
    auto p = stage3_params();
    p.mesh.INDEX_Serial = 0;
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 2);
    EXPECT_THROW((Transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "")), std::invalid_argument);
}

TEST(ModuleInterfaces, WisdomExportImportAndMissingFileContract) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    const auto file = directory.path / "plans.wisdom";
    EXPECT_FALSE(wisdom::import(file.string()));
    { Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, file.string()); }
    EXPECT_TRUE(std::filesystem::is_regular_file(file));
    EXPECT_FALSE(std::filesystem::exists(file.string()+".tmp"));
    EXPECT_TRUE(wisdom::import(file.string()));
    test_support::write(file, "invalid FFTW wisdom");
    EXPECT_FALSE(wisdom::import(file.string()));
    EXPECT_NO_THROW(wisdom::export_wisdom((directory.path / "missing/file").string()));
}

TEST(ModuleInterfaces, VectorOperationsRejectInvalidAxesAndCompareComplexFields) {
    const auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms t(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    auto a = make_scalar_field(m), b = make_scalar_field(m);
    for (int axis : {-1, 3}) EXPECT_THROW(vector_ops::derivative(axis, b, a, m, t), std::invalid_argument);
    a.at_z(-m.NZS, 0, 0, 0) = {3.,4.};
    EXPECT_DOUBLE_EQ(vector_ops::max_distance(a, b, m), 5.);
    EXPECT_DOUBLE_EQ(vector_ops::max_distance(a, a, m), 0.);
}
