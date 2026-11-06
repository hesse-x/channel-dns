#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

#include "core/params.hpp"
#include "test_support.hpp"

using test_support::replace;

class TomlConfig : public ::testing::Test {
protected:
    test_support::TempDirectory directory;
    std::string sample = test_support::read(CHANNEL_DNS_SOURCE_DIR "/examples/config.toml");
    Params load(const std::string& text) {
        test_support::write(directory.path / "config.toml", text);
        return Params::load(directory.path.string());
    }
};

TEST_F(TomlConfig, SupportsCommentsMultilineArraysIntegerRealsAndNumericSeparators) {
    auto text = replace(sample, "Re = 2.8e3", "Re = 2_800 # integer accepted as real");
    text = replace(text, "a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = [\n  1, # first\n  2.5,\n  -3e-2,\n]");
    text = replace(text, "b_V = [0.0e0, 0.0e0, 0.0e0]", "b_V = [4, 5, 6]");
    const auto p = load(text);
    EXPECT_DOUBLE_EQ(p.parm.Re, 2800);
    EXPECT_EQ(p.parm.a_V, (std::array<double, 3>{1, 2.5, -.03}));
    EXPECT_EQ(p.parm.b_V, (std::array<double, 3>{4, 5, 6}));
}

TEST_F(TomlConfig, OutputModesAndOptionalOutputTable) {
    for (const auto& [name, mode] : {std::pair{"ascii", OutputFormat::Ascii},
                                  {"binary", OutputFormat::Binary}, {"both", OutputFormat::Both}}) {
        SCOPED_TRACE(name);
        EXPECT_EQ(load(replace(sample, "format = \"ascii\"", std::string("format = \"") + name + "\"")).output_format, mode);
    }
    EXPECT_EQ(load(sample.substr(0, sample.find("[output]"))).output_format, OutputFormat::Ascii);
}

TEST_F(TomlConfig, IndependentOutputIntervalsAreOptionalAndValidated) {
    const auto defaults = load(sample);
    EXPECT_EQ(defaults.physical_interval, -1);
    EXPECT_EQ(defaults.statistics_interval, -1);
    const auto p = load(sample + "\nphysical_interval = 0\nstatistics_interval = 10000\n");
    EXPECT_EQ(p.physical_interval, 0);
    EXPECT_EQ(p.statistics_interval, 10000);
    for (const auto bad : {"physical_interval = -1", "physical_interval = 1.5",
                           "statistics_interval = 0", "statistics_interval = -1"})
        EXPECT_THROW(load(sample + "\n" + bad + "\n"), std::runtime_error);
}

TEST_F(TomlConfig, AcceptsCaseInsensitiveNamesAndOutputModes) {
    for (const bool upper : {false, true}) {
        auto text = sample;
        for (char& c : text) {
            if (upper && c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
            if (!upper && c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
        }
        const auto p = load(text);
        EXPECT_EQ(p.mesh.NX, 32);
        EXPECT_DOUBLE_EQ(p.parm.Re, 2800);
        EXPECT_EQ(p.parm.a_V, (std::array<double, 3>{1, 1, 1}));
        EXPECT_EQ(p.output_format, OutputFormat::Ascii);
    }
    auto text = replace(sample, "[mesh]", "[MeSh]");
    text = replace(text, "INDEX_Serial", "InDeX_sErIaL");
    text = replace(text, "format = \"ascii\"", "FoRmAt = \"BiNaRy\"");
    EXPECT_EQ(load(text).output_format, OutputFormat::Binary);
    EXPECT_EQ(load(replace(text, "BiNaRy", "BoTh")).output_format, OutputFormat::Both);
}

TEST_F(TomlConfig, DerivedQuantitiesCoverBothPressureModesAndNonzeroAngle) {
    auto text = replace(sample, "THETA = 0.0e0", "THETA = 60");
    auto p = load(text);
    EXPECT_NEAR(p.THETA_RAD, std::acos(-1.0) / 3, 1e-15);
    EXPECT_NEAR(p.QXmean, .5, 1e-15);
    EXPECT_NEAR(p.QZmean, std::sqrt(3.) / 2, 1e-15);
    EXPECT_EQ(p.GRAD_PRE_X, 0);
    EXPECT_EQ(p.GRAD_PRE_Z, 0);
    p = load(replace(text, "INDEX_GRAD_PRE = 0", "INDEX_GRAD_PRE = 1"));
    EXPECT_NEAR(p.GRAD_PRE_X, -.5, 1e-15);
    EXPECT_NEAR(p.GRAD_PRE_Z, -std::sqrt(3.) / 2, 1e-15);
    EXPECT_EQ(p.QXmean, 0);
    EXPECT_EQ(p.QZmean, 0);
}

TEST_F(TomlConfig, TomlTakesPrecedenceAndInvalidTomlDoesNotFallBack) {
    for (const auto name : {"meshparm", "parm", "control"})
        test_support::write(directory.path / (std::string(name) + "_channel_with_Microblow.txt"), "invalid text configuration");
    EXPECT_EQ(load(sample).mesh.NX, 32);
    EXPECT_THROW(load("[mesh\n"), std::exception);
}

TEST_F(TomlConfig, EveryRequiredKeyMustBePresent) {
    // Derive deletions from the public example so additions become required coverage.
    std::size_t start = 0;
    while (start < sample.find("[output]")) {
        const auto end = sample.find('\n', start);
        const auto line = sample.substr(start, end - start);
        if (line.find('=') != std::string::npos) {
            SCOPED_TRACE(line);
            auto missing = sample;
            missing.erase(start, end - start + 1);
            EXPECT_THROW(load(missing), std::exception);
        }
        start = end + 1;
    }
}

class InvalidTomlValue : public TomlConfig, public ::testing::WithParamInterface<std::pair<std::string, std::string>> {};
TEST_P(InvalidTomlValue, RejectsInvalidInput) {
    const auto& [before, after] = GetParam();
    EXPECT_THROW(load(replace(sample, before, after)), std::exception);
}
INSTANTIATE_TEST_SUITE_P(Contract, InvalidTomlValue, ::testing::ValuesIn(std::vector<std::pair<std::string,std::string>>{
    {"NX = 32", "NX = true"},
    {"NX = 32", "NX = \"32\""}, {"NX = 32", "NX = 32.0"},
    {"NX = 32", "NX = 2147483648"}, {"NX = 32", "NX = -2147483649"},
    {"NX = 32", "NX = 9223372036854775808"},
    {"DT = 3.0e-3", "DT = \"0.003\""}, {"DT = 3.0e-3", "DT = true"},
    {"DT = 3.0e-3", "DT = nan"}, {"DT = 3.0e-3", "DT = inf"},
    {"DT = 3.0e-3", "DT = -inf"}, {"DT = 3.0e-3", "DT = 3.0D-3"},
    {"a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = []"},
    {"a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = [1,2]"},
    {"a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = [1,2,3,4]"},
    {"a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = [1,true,3]"},
    {"a_V = [1.0e0, 1.0e0, 1.0e0]", "a_V = [1,nan,3]"},
    {"b_V = [0.0e0, 0.0e0, 0.0e0]", "b_V = [1,inf,3]"},
    {"b_V = [0.0e0, 0.0e0, 0.0e0]", "b_V = 1"},
    {"format = \"ascii\"", "format = \"typo\""},
    {"format = \"ascii\"", "format = true"},
    {"format = \"ascii\"", ""},
    {"NX = 32", "NX = 32\nNX = 32"},
    {"NX = 32", "NX = 32\nnx = 64"},
    {"[mesh]", "[MESH]\n[mesh]"},
    {"NX = 32", "NX = 32\nunknown = 0"},
    {"Re = 2.8e3", "Re = 2.8e3\nunknown = 0"},
    {"Nstat = 50000", "Nstat = 50000\nunknown = 0"},
    {"format = \"ascii\"", "format = \"ascii\"\nunknown = 0"},
    {"[mesh]", "[meshh]"}, {"[parm]", "[parameters]"},
    {"[control]", "[controls]"}, {"[output]", "[outputs]"}
}));
