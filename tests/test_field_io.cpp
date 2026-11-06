#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>

#include "core/mesh.hpp"
#include "fft/transforms.hpp"
#include "io/field_io.hpp"
#include "stage3_fixture.hpp"

namespace {
std::filesystem::path temp_path() {
    return std::filesystem::temp_directory_path() /
        ("channel-dns-io-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
}
}

TEST(FieldIO, RestartRoundTripIsBitExact) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    const auto dir = temp_path();
    std::filesystem::create_directories(dir);
    FieldIO io(m, transforms, dir);
    VectorField source = make_vector_field(m), restored = make_vector_field(m);
    for (int i = 0; i < source.size(); ++i) source.data()[i] = {i * 0.125, -i * 0.25};
    io.write_restart(17, 0.75, source);
    std::int64_t nt = 0;
    double time = 0.0;
    io.read_restart(nt, time, restored);
    EXPECT_EQ(nt, 17);
    EXPECT_EQ(time, 0.75);
    for (int i = 0; i < source.size(); ++i) EXPECT_EQ(restored.data()[i], source.data()[i]);
    EXPECT_FALSE(std::filesystem::exists(dir / "field_S_ 0.dat.tmp"));
    std::filesystem::remove_all(dir);
}

TEST(FieldIO, EktAppends) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    const auto dir = temp_path();
    std::filesystem::create_directories(dir);
    FieldIO io(m, transforms, dir);
    PhysicalSnapshot snapshot;
    for (auto& v : snapshot.velocity) v.assign(4, 0.0);
    snapshot.velocity[1].assign(4, 2.0);
    io.append_ekt(1, 0.1, snapshot);
    io.append_ekt(2, 0.2, snapshot);
    std::ifstream in(dir / "EKT 0.dat");
    double t1, e1, t2, e2;
    ASSERT_TRUE(in >> t1 >> e1 >> t2 >> e2);
    EXPECT_DOUBLE_EQ(e1, 16.0);
    EXPECT_DOUBLE_EQ(e2, 16.0);
    std::filesystem::remove_all(dir);
}
TEST(FieldIO, BinaryRestartRoundTripIsBitExact) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    const auto dir = temp_path();
    std::filesystem::create_directories(dir);
    FieldIO io(m, transforms, dir, nullptr, OutputFormat::Binary);
    VectorField source = make_vector_field(m), restored = make_vector_field(m);
    for (int i = 0; i < source.size(); ++i) source.data()[i] = {i * 0.125, -i * 0.25};
    io.write_restart(17, 0.75, source);
    std::int64_t nt = 0;
    double time = 0.0;
    io.read_restart(nt, time, restored);
    EXPECT_EQ(nt, 17);
    EXPECT_EQ(time, 0.75);
    for (int i = 0; i < source.size(); ++i) EXPECT_EQ(restored.data()[i], source.data()[i]);
    EXPECT_FALSE(std::filesystem::exists(dir / "field_S_ 0.dat.tmp"));
    std::filesystem::remove_all(dir);
}


TEST(FieldIO, PhysicalTextUsesComponentFirstOrder) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    const auto dir = temp_path();
    std::filesystem::create_directories(dir);
    FieldIO io(m, transforms, dir);
    PhysicalSnapshot snapshot;
    for (int c = 0; c < 3; ++c) snapshot.velocity[c] = {10.0 + c, 20.0 + c};
    snapshot.pressure = {30.0, 31.0};
    io.write_physical(9, 0.5, snapshot);
    std::ifstream in(dir / "field_P_ 0.dat");
    double v;
    for (double expected : {9., .5, 10., 11., 12., 20., 21., 22., 30., 31.}) {
        ASSERT_TRUE(in >> v);
        EXPECT_EQ(v, expected);
    }
    std::filesystem::remove_all(dir);
}

TEST(FieldIO, TruncatedRestartLeavesStateUnchanged) {
    Params p = stage3_params();
    Mesh m(p.mesh, p.parm.XL, p.parm.ZL, 0, 1);
    Transforms transforms(m, p.parm.XL, p.parm.ZL, p.parm.Re, "");
    const auto dir = temp_path();
    std::filesystem::create_directories(dir);
    auto velocity = make_vector_field(m);
    velocity.fill({3.0, 4.0});
    for (auto format : {OutputFormat::Ascii, OutputFormat::Binary}) {
        FieldIO io(m, transforms, dir, nullptr, format);
        io.write_restart(17, 0.75, velocity);
        const auto path = dir / (format == OutputFormat::Binary ? "field_S_ 0.bin" : "field_S_ 0.dat");
        std::filesystem::resize_file(path, 64);
        std::int64_t nt = 99;
        double time = 1.0;
        EXPECT_THROW(io.read_restart(nt, time, velocity), std::runtime_error);
        EXPECT_EQ(nt, 99);
        EXPECT_EQ(time, 1.0);
        for (int i = 0; i < velocity.size(); ++i) EXPECT_EQ(velocity.data()[i], std::complex<double>(3.,4.));
    }
    std::filesystem::remove_all(dir);
}

#include <limits>
#include "test_support.hpp"

class FieldIOContract : public ::testing::Test {
protected:
    test_support::TempDirectory directory;
    Params p = stage3_params();
    Mesh mesh{p.mesh, p.parm.XL, p.parm.ZL, 0, 1};
    Transforms transforms{mesh, p.parm.XL, p.parm.ZL, p.parm.Re, ""};
    VectorField source = make_vector_field(mesh);
    void SetUp() override {
        for (int i = 0; i < source.size(); ++i) source.data()[i] = {i*.125, -.25*i};
    }
};

TEST_F(FieldIOContract, SupportsBothModesFallbackAndOldBinaryDatFiles) {
    FieldIO writer(mesh, transforms, directory.path, nullptr, OutputFormat::Both);
    writer.write_restart(7, .25, source);
    EXPECT_TRUE(std::filesystem::exists(directory.path / "field_S_ 0.dat"));
    EXPECT_TRUE(std::filesystem::exists(directory.path / "field_S_ 0.bin"));
    for (auto format : {OutputFormat::Ascii, OutputFormat::Binary, OutputFormat::Both}) {
        FieldIO reader(mesh, transforms, directory.path, nullptr, format);
        auto restored = make_vector_field(mesh);
        std::int64_t nt = 0; double time = 0;
        reader.read_restart(nt, time, restored);
        EXPECT_EQ(nt, 7);
        EXPECT_EQ(time, .25);
        for (int i = 0; i < source.size(); ++i) EXPECT_EQ(restored.data()[i], source.data()[i]);
    }
    std::filesystem::remove(directory.path / "field_S_ 0.dat");
    FieldIO text_reader(mesh, transforms, directory.path);
    auto restored = make_vector_field(mesh);
    std::int64_t nt = 0; double time = 0;
    EXPECT_NO_THROW(text_reader.read_restart(nt, time, restored));
    std::filesystem::rename(directory.path / "field_S_ 0.bin", directory.path / "field_S_ 0.dat");
    FieldIO binary_reader(mesh, transforms, directory.path, nullptr, OutputFormat::Binary);
    EXPECT_NO_THROW(binary_reader.read_restart(nt, time, restored));
    EXPECT_EQ(nt, 7);
    for (int i = 0; i < source.size(); ++i) EXPECT_EQ(restored.data()[i], source.data()[i]);
}

TEST_F(FieldIOContract, ToPhysicalPreservesInputsAndSeparatesComponents) {
    FieldIO io(mesh, transforms, directory.path);
    source.fill({});
    for (int c = 0; c < 3; ++c) source.at_z(0, 0, 0, c) = 2.+c;
    auto pressure = make_scalar_field(mesh);
    pressure.at_z(0, 0, 0, 0) = 7.;
    const auto original = source;
    const auto snapshot = io.to_physical(source, pressure);
    const std::size_t expected_size = mesh.NXTNZ*(mesh.NYP+1);
    EXPECT_EQ(snapshot.pressure.size(), expected_size);
    for (int c = 0; c < 3; ++c) {
        ASSERT_EQ(snapshot.velocity[c].size(), expected_size);
        for (const double value : snapshot.velocity[c]) EXPECT_NEAR(value, 2.+c, 1e-13);
    }
    for (const double value : snapshot.pressure) EXPECT_NEAR(value, 7., 1e-13);
    for (int i = 0; i < source.size(); ++i) EXPECT_EQ(source.data()[i], original.data()[i]);
    EXPECT_EQ(pressure.at_z(0, 0, 0, 0), std::complex<double>(7.));
    source.fill({}); pressure.fill({});
    const auto& zero = io.to_physical(source, pressure);
    for (const auto& component : zero.velocity) for (const double value : component) EXPECT_EQ(value, 0.);
    for (const double value : zero.pressure) EXPECT_EQ(value, 0.);
}

TEST_F(FieldIOContract, OutputFailureLeavesNoTemporaryFiles) {
    // A missing parent makes failure deterministic even when tests run as root.
    FieldIO io(mesh, transforms, directory.path / "missing", nullptr, OutputFormat::Both);
    EXPECT_THROW(io.write_restart(1, .1, source), std::runtime_error);
    PhysicalSnapshot snapshot;
    EXPECT_THROW(io.write_physical(1, .1, snapshot), std::runtime_error);
    EXPECT_THROW(io.append_ekt(1, .1, snapshot), std::runtime_error);
    EXPECT_TRUE(std::filesystem::is_empty(directory.path));
}

TEST_F(FieldIOContract, TextParserHandlesDExponentsCommasAndRejectsCorruptPayloads) {
    FieldIO io(mesh, transforms, directory.path);
    io.write_restart(7, .25, source);
    const auto path = directory.path / "field_S_ 0.dat";
    const auto original = test_support::read(path);
    auto text_variant = original;
    for (char& c : text_variant) if (c == 'e') c = 'D';
    text_variant = test_support::replace(text_variant, "7 ", "7, ");
    test_support::write(path, text_variant);
    auto restored = make_vector_field(mesh);
    std::int64_t nt = 99; double time = 99;
    ASSERT_NO_THROW(io.read_restart(nt, time, restored));
    for (int i = 0; i < source.size(); ++i) EXPECT_EQ(source.data()[i], restored.data()[i]);
    for (const auto& corrupt : {original + "1\n", test_support::replace(original, "(", "["),
                              test_support::replace(original, ",", ";"), std::string("-1 0\n"),
                              test_support::replace(original, "(", "(nan,")}) {
        test_support::write(path, corrupt);
        restored.fill({9.,8.}); nt = 99; time = 99;
        EXPECT_THROW(io.read_restart(nt, time, restored), std::exception);
        EXPECT_EQ(nt, 99);
        EXPECT_EQ(time, 99);
        for (int i = 0; i < restored.size(); ++i) EXPECT_EQ(restored.data()[i], std::complex<double>(9.,8.));
    }
}

class CorruptBinaryHeader : public FieldIOContract, public ::testing::WithParamInterface<std::size_t> {};
TEST_P(CorruptBinaryHeader, RejectsIncompatibleMetadataWithoutChangingState) {
    FieldIO io(mesh, transforms, directory.path, nullptr, OutputFormat::Binary);
    io.write_restart(7, .25, source);
    const auto path = directory.path / "field_S_ 0.bin";
    auto bytes = test_support::read(path);
    bytes[GetParam()] ^= 0x40;
    test_support::write(path, bytes);
    auto restored = make_vector_field(mesh);
    restored.fill({9.,8.});
    std::int64_t nt = 99; double time = 99;
    EXPECT_THROW(io.read_restart(nt, time, restored), std::exception);
    EXPECT_EQ(nt, 99);
    EXPECT_EQ(time, 99);
    for (int i = 0; i < restored.size(); ++i) EXPECT_EQ(restored.data()[i], std::complex<double>(9.,8.));
}
INSTANTIATE_TEST_SUITE_P(HeaderFields, CorruptBinaryHeader,
    ::testing::Values(std::size_t{0}, 8, 12, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120));
