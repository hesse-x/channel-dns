#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>

#include "core/mesh.hpp"
#include "sim/statistics.hpp"
#include "stage3_fixture.hpp"

TEST(Statistics, UpdatesMeansWithoutMutatingPressure) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Statistics statistics(m);
    const std::size_t count = static_cast<std::size_t>(m.NZ) * (m.NYP + 1) * m.NX;
    PhysicalVectorField velocity;
    velocity[0].assign(count, 1.0);
    velocity[1].assign(count, 2.0);
    velocity[2].assign(count, 3.0);
    PhysicalScalarField pressure(count, 10.0);
    const auto original = pressure;
    statistics.accumulate(velocity, pressure, -1.0, 2.0);
    statistics.accumulate(velocity, pressure, 3.0, 4.0);
    EXPECT_EQ(statistics.count(), 2);
    EXPECT_DOUBLE_EQ(statistics.dpdx_mean(), 1.0);
    EXPECT_DOUBLE_EQ(statistics.dpdz_mean(), 3.0);
    EXPECT_EQ(pressure, original);

    const auto directory = std::filesystem::temp_directory_path() /
        ("channel-dns-stat-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    statistics.write(directory);
    EXPECT_TRUE(std::filesystem::is_regular_file(directory / "Result/0/DPDX_DPDZ.dat"));
    EXPECT_TRUE(std::filesystem::is_regular_file(directory / "Result/4/V_Pre_moment4_Field_ 0.dat"));
    EXPECT_FALSE(std::filesystem::exists(directory / "Result/0/DPDX_DPDZ.dat.tmp"));
    std::filesystem::remove_all(directory);
}

#include <fstream>
#include <sstream>
#include "test_support.hpp"

TEST(Statistics, RejectsWrongFieldSizesWithoutChangingCountOrMeans) {
    auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Statistics stats(m);
    const std::size_t n = m.NXTNZ * (m.NYP+1);
    PhysicalVectorField velocity;
    for (auto& c : velocity) c.assign(n, 1.);
    PhysicalScalarField pressure(n, 2.);
    stats.accumulate(velocity, pressure, 3., 4.);
    for (int component = 0; component < 4; ++component) {
        auto bad_v = velocity;
        auto bad_p = pressure;
        if (component == 3) bad_p.pop_back(); else bad_v[component].pop_back();
        EXPECT_THROW(stats.accumulate(bad_v, bad_p, 99., 99.), std::runtime_error);
        EXPECT_EQ(stats.count(), 1);
        EXPECT_EQ(stats.dpdx_mean(), 3.);
        EXPECT_EQ(stats.dpdz_mean(), 4.);
    }
}

TEST(Statistics, WritesAnalyticMeansCrossMomentsAndPressureCorrection) {
    test_support::TempDirectory directory;
    auto p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Statistics stats(m);
    const std::size_t n = m.NXTNZ * (m.NYP+1);
    PhysicalVectorField v;
    PhysicalScalarField pressure(n);
    // Two uniform samples: v=(1,2,3), p_corrected=4 and v=(3,4,5), p_corrected=6.
    for (int sample = 0; sample < 2; ++sample) {
        for (int c = 0; c < 3; ++c) v[c].assign(n, 1+c+2*sample);
        const double kinetic = .5 * (v[0][0]*v[0][0]+v[1][0]*v[1][0]+v[2][0]*v[2][0]);
        pressure.assign(n, 4+2*sample+kinetic);
        stats.accumulate(v, pressure, -1+2*sample, 2+2*sample);
    }
    stats.write(directory.path);
    const std::vector<std::vector<double>> expected{
        {2,3,4,5}, {5,10,17,7,13,9,26,11,16,21}, {14,36,76}, {41,136,353}};
    const std::vector<std::string> names{"mean", "moment2", "moment3", "moment4"};
    for (int moment = 0; moment < 4; ++moment) {
        const auto path = directory.path / "Result" / std::to_string(moment+1) /
            ("V_Pre_"+names[moment]+"_Field_ 0.dat");
        std::ifstream in(path);
        ASSERT_TRUE(in);
        std::string line;
        for (int iy = 0; iy <= m.NYP; ++iy) {
            ASSERT_TRUE(std::getline(in, line));
            EXPECT_EQ(line.size(), 19*(expected[moment].size()+1));
            std::istringstream row(line);
            double value;
            ASSERT_TRUE(row >> value);
            EXPECT_NEAR(value, m.YC[iy], 5.1e-9);
            for (const auto wanted : expected[moment]) {
                ASSERT_TRUE(row >> value);
                EXPECT_DOUBLE_EQ(value, wanted);
            }
            EXPECT_FALSE(row >> value);
        }
        EXPECT_FALSE(std::getline(in, line));
    }
    std::ifstream gradients(directory.path / "Result/0/DPDX_DPDZ.dat");
    int count; double x, z;
    ASSERT_TRUE(gradients >> count >> x >> z);
    EXPECT_EQ(count, 2);
    EXPECT_EQ(x, 0);
    EXPECT_EQ(z, 3);
}
