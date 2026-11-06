// Params 单测：内嵌参数文件样例，覆盖 D 指数、逗号/空格混合分隔、
// int 遇小数点报错、EOF/缺文件报错。
#include "core/params.hpp"

#include <cmath>
#include <filesystem>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

namespace {

std::string make_dir() {
    std::string dir = ::testing::TempDir() + "params_test";
    std::string cmd = "mkdir -p " + dir;
    if (std::system(cmd.c_str()) != 0) ADD_FAILURE() << "mkdir failed";
    return dir;
}

void write(const std::string& path, const std::string& content) {
    std::ofstream o(path);
    o << content;
}

// 与真实文件同构的样例：D 指数、逗号/空格混合、无空格逗号分隔
void write_standard_set(const std::string& d) {
    write(d + "/meshparm_channel_with_Microblow.txt",
          "NX,NY,NZ,INDEX_disturb_Z\n"
          "32,32,32,3\n"
          "NXS,NYS,NZS\n"
          "16, 30, 16\n"
          "NX2,NY2,NZ2\n"
          "32 32 32\n"
          "INDEX_Serial\n"
          "1\n");
    write(d + "/parm_channel_with_Microblow.txt",
          "XL,ZL,THETA\n"
          "1.0D1,1.0D1,4.5D1\n"
          "INDEX_GRAD_PRE,Re\n"
          "0,2.8D3\n"
          "INDEX_UNSTEADY_BC\n"
          "0\n"
          "G_pre,a_pre,b_pre\n"
          "1.0D0 3.0D-3,2.0D0\n"
          "a_V(I),b_V(I),I=1,3\n"
          "1.0D0,0.0D0\n"
          "5.0D-1, 1.0D0\n"
          "2.5D0,3.0D0\n");
    write(d + "/control_channel_with_Microblow.txt",
          "INDEX_CONT,EPS_init\n"
          "0,0.1D0\n"
          "NSTEPS,NDUMFRQ,NPRNFRQ,DT\n"
          "5,200000,20000,3.0D-3\n"
          "NCOUNT_disturb,NDISTURBFRQ,EPS_disturb\n"
          "30,300,0.1D0\n"
          "V_Microblow0\n"
          "1.0D-3\n"
          "Nstat\n"
          "50000\n");
}

} // namespace

TEST(Params, ParsesStandardFiles) {
    const std::string d = make_dir() + "_std";
    std::system(("mkdir -p " + d).c_str());
    write_standard_set(d);
    Params p = Params::load(d);

    EXPECT_EQ(p.mesh.NX, 32);
    EXPECT_EQ(p.mesh.INDEX_disturb_Z, 3);
    EXPECT_EQ(p.mesh.NYS, 30);
    EXPECT_EQ(p.mesh.NX2, 32);
    EXPECT_EQ(p.mesh.INDEX_Serial, 1);

    EXPECT_DOUBLE_EQ(p.parm.XL, 10.0);
    EXPECT_DOUBLE_EQ(p.parm.THETA, 45.0);
    EXPECT_EQ(p.parm.INDEX_GRAD_PRE, 0);
    EXPECT_DOUBLE_EQ(p.parm.Re, 2800.0);
    EXPECT_DOUBLE_EQ(p.parm.G_pre, 1.0);
    EXPECT_DOUBLE_EQ(p.parm.a_pre0, 3.0e-3);
    EXPECT_DOUBLE_EQ(p.parm.b_pre0, 2.0);
    EXPECT_DOUBLE_EQ(p.parm.a_V[1], 0.5);   // 逗号+空格混合分隔
    EXPECT_DOUBLE_EQ(p.parm.b_V[2], 3.0);

    EXPECT_EQ(p.ctrl.NSTEPS, 5);
    EXPECT_DOUBLE_EQ(p.ctrl.DT, 3.0e-3);
    EXPECT_EQ(p.ctrl.NCOUNT_disturb, 30);
    EXPECT_DOUBLE_EQ(p.ctrl.V_Microblow0, 1.0e-3);
    EXPECT_EQ(p.ctrl.Nstat, 50000);

    // 派生量：INDEX_GRAD_PRE==0 → QXmean/QZmean 分支
    const double pi = std::acos(-1.0);
    EXPECT_DOUBLE_EQ(p.THETA_RAD, 45.0 * pi / 180.0);
    EXPECT_DOUBLE_EQ(p.QXmean, std::cos(p.THETA_RAD));
    EXPECT_DOUBLE_EQ(p.QZmean, std::sin(p.THETA_RAD));
}

TEST(Params, GradPreBranch) {
    const std::string d = make_dir() + "_gp";
    std::system(("mkdir -p " + d).c_str());
    write_standard_set(d);
    // 改成 INDEX_GRAD_PRE=1
    write(d + "/parm_channel_with_Microblow.txt",
          "XL,ZL,THETA\n"
          "1.0D1,1.0D1,9.0D1\n"
          "INDEX_GRAD_PRE,Re\n"
          "1,2.8D3\n"
          "INDEX_UNSTEADY_BC\n"
          "0\n"
          "G_pre,a_pre,b_pre\n"
          "0.0D0,0.0D0,1.0D0\n"
          "a_V(I),b_V(I),I=1,3\n"
          "1.0D0,0.0D0\n1.0D0,0.0D0\n1.0D0,0.0D0\n");
    Params p = Params::load(d);
    const double pi = std::acos(-1.0);
    EXPECT_DOUBLE_EQ(p.GRAD_PRE_X, -std::cos(pi / 2.0));
    EXPECT_DOUBLE_EQ(p.GRAD_PRE_Z, -std::sin(pi / 2.0));
}

TEST(Params, IntWithDecimalPointThrows) {
    const std::string d = make_dir() + "_badint";
    std::system(("mkdir -p " + d).c_str());
    write_standard_set(d);
    write(d + "/meshparm_channel_with_Microblow.txt",
          "NX,NY,NZ,INDEX_disturb_Z\n"
          "32.0,32,32,3\n"           // int 位置出现小数点 → 必须报错
          "NXS,NYS,NZS\n16,30,16\n"
          "NX2,NY2,NZ2\n32,32,32\n"
          "INDEX_Serial\n1\n");
    EXPECT_THROW(Params::load(d), std::runtime_error);
}

TEST(Params, MissingFileThrows) {
    EXPECT_THROW(Params::load("/nonexistent_dir_xyz"), std::runtime_error);
}

TEST(Params, TokenCountMismatchThrows) {
    const std::string d = make_dir() + "_badcnt";
    std::system(("mkdir -p " + d).c_str());
    write_standard_set(d);
    write(d + "/control_channel_with_Microblow.txt",
          "INDEX_CONT,EPS_init\n"
          "0,0.1D0,9.9D9\n"          // 多一个 token → 报错
          "NSTEPS,NDUMFRQ,NPRNFRQ,DT\n5,200000,20000,3.0D-3\n"
          "NCOUNT_disturb,NDISTURBFRQ,EPS_disturb\n30,300,0.1D0\n"
          "V_Microblow0\n1.0D-3\n"
          "Nstat\n50000\n");
    EXPECT_THROW(Params::load(d), std::runtime_error);
}

TEST(Params, ExampleTomlParameters) {
    const auto p = Params::load(std::string(CHANNEL_DNS_SOURCE_DIR) + "/examples");
    EXPECT_EQ(p.mesh.NX, 32);
    EXPECT_EQ(p.mesh.NY, 32);
    EXPECT_EQ(p.mesh.NZ, 32);
    EXPECT_EQ(p.mesh.NXS, 16);
    EXPECT_EQ(p.mesh.NYS, 30);
    EXPECT_EQ(p.mesh.NZS, 16);
    EXPECT_EQ(p.mesh.INDEX_Serial, 1);
    EXPECT_DOUBLE_EQ(p.parm.XL, 10.0);
    EXPECT_DOUBLE_EQ(p.parm.ZL, 10.0);
    EXPECT_DOUBLE_EQ(p.parm.Re, 2800.0);
    EXPECT_EQ(p.ctrl.NSTEPS, 5);
    EXPECT_DOUBLE_EQ(p.ctrl.DT, 0.003);
    EXPECT_EQ(p.output_format, OutputFormat::Ascii);
}

TEST(Params, TomlRejectsMissingAndUnknownKeys) {
    const auto dir = make_dir() + "_toml";
    std::filesystem::create_directories(dir);
    std::ifstream sample(std::string(CHANNEL_DNS_SOURCE_DIR) + "/examples/config.toml");
    const std::string original((std::istreambuf_iterator<char>(sample)), {});
    write(dir + "/config.toml", original + "typo = 1\n");
    EXPECT_THROW(Params::load(dir), std::runtime_error);
    auto missing = original;
    const auto start = missing.find("NX =");
    missing.erase(start, missing.find('\n', start) - start + 1);
    write(dir + "/config.toml", missing);
    EXPECT_THROW(Params::load(dir), std::runtime_error);
    write(dir + "/config.toml", original.substr(0, original.find("format =")) + "format = \"both\"\n");
    EXPECT_EQ(Params::load(dir).output_format, OutputFormat::Both);
    auto fractional = original;
    fractional.replace(fractional.find("NX = 32"), 7, "NX = 32.0");
    write(dir + "/config.toml", fractional);
    EXPECT_THROW(Params::load(dir), std::runtime_error);
    auto wrong_array = original;
    const auto array_begin = wrong_array.find("a_V =");
    wrong_array.replace(array_begin, wrong_array.find('\n', array_begin) - array_begin, "a_V = [1.0, 2.0]");
    write(dir + "/config.toml", wrong_array);
    EXPECT_THROW(Params::load(dir), std::runtime_error);
}
