// 3D 装配测试（step2 §7.2 F1-F5）：经 Transforms 全装配验证 xzfft / channel_xyzfft /
// channel_xyzfft2 串行路径。两配置：
//   C-FULL  （NXS=NZS=16=NXH/NZH，NYS=30）—— 满谱 Nyquist 分支（生产配置，LIVE 代码）
//   C-TRUNC （NXS=NZS=10，NYS=20）—— 截断分支
// 本配置粗/细同尺寸（NX2=NY2=NZ2=32=NX=NY=NZ），故三组变换同布局，分别自洽验证。
//
// 期望值构造独立于 FFTW：x/z 用 cos 单模态手算（r2c/c2c 半峰），y 用 DCT-I 定义式
// （dct1_definition）按生产归一化（[0]/(2NY)、[1..NYS]/NY）全谱计算 —— 不假设峰型。
#include <gtest/gtest.h>
#include <fftw3.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include "core/mesh.hpp"
#include "core/params.hpp"
#include "core/cpu_device.hpp"
#include "fft/transforms.hpp"

using cdouble = std::complex<double>;

namespace {
constexpr double PI = 3.14159265358979323846;

MeshParm make_meshparm(int NXS, int NYS, int NZS) {
    MeshParm mp{};
    mp.NX = 32; mp.NY = 32; mp.NZ = 32;
    mp.NXS = NXS; mp.NYS = NYS; mp.NZS = NZS;
    mp.NX2 = 32; mp.NY2 = 32; mp.NZ2 = 32;
    mp.INDEX_disturb_Z = 0;
    mp.INDEX_Serial = 1;
    return mp;
}

// 串行下 U 谱场布局 [iz+NZS][iy][kx]（NYS+1 × NXSP+1），iz ∈ -NZS..NZS。
struct Spectral {
    const Mesh& m;
    std::vector<cdouble> v;
    explicit Spectral(const Mesh& mesh) : m(mesh),
        v(static_cast<size_t>(mesh.NZS2P1) * (mesh.NYS + 1) * (mesh.NXSP + 1), {0, 0}) {}
    size_t idx(int iz, int iy, int kx) const {
        return (static_cast<size_t>(iz + m.NZS) * (m.NYS + 1) + iy) * (m.NXSP + 1) + kx;
    }
    cdouble& at(int iz, int iy, int kx) { return v[idx(iz, iy, kx)]; }
    cdouble  at(int iz, int iy, int kx) const { return v[idx(iz, iy, kx)]; }
};

// 物理场 [iz][iy][ix]（串行 NYP=NY，NYP+1=NY+1）。
struct Physical {
    const Mesh& m;
    std::vector<double> v;
    explicit Physical(const Mesh& mesh) : m(mesh),
        v(static_cast<size_t>(mesh.NZ) * (mesh.NYP + 1) * mesh.NX, 0.0) {}
    size_t idx(int iz, int iy, int ix) const {
        return (static_cast<size_t>(iz) * (m.NYP + 1) + iy) * m.NX + ix;
    }
    double& at(int iz, int iy, int ix) { return v[idx(iz, iy, ix)]; }
    double  at(int iz, int iy, int ix) const { return v[idx(iz, iy, ix)]; }
};

double maxdiff(const Spectral& a, const Spectral& b) {
    double md = 0.0;
    for (size_t i = 0; i < a.v.size(); ++i) md = std::max(md, std::abs(a.v[i] - b.v[i]));
    return md;
}
double maxdiff(const Physical& a, const Physical& b) {
    double md = 0.0;
    for (size_t i = 0; i < a.v.size(); ++i) md = std::max(md, std::fabs(a.v[i] - b.v[i]));
    return md;
}

// DCT-I（REDFT00）定义式（§0.1-5）：
//   Y_k = x_0 + (-1)^k·x_{n-1} + 2·Σ_{j=1}^{n-2} x_j·cos(πjk/(n-1))
std::vector<double> dct1_definition(const std::vector<double>& x) {
    const int n = static_cast<int>(x.size());
    const int ny = n - 1;
    std::vector<double> y(n, 0.0);
    for (int k = 0; k < n; ++k) {
        double s = x[0] + ((k % 2 == 0) ? x[ny] : -x[ny]);
        for (int j = 1; j <= ny - 1; ++j) s += 2.0 * x[j] * std::cos(PI * j * k / ny);
        y[k] = s;
    }
    return y;
}

// 物理 y 线 T_my 经生产 'F' 的期望 y 谱：DCT-I 定义式 → [0]/(2NY)、[1..NYS]/NY，其余 0。
std::vector<double> y_spectrum_expected(int my, int NY, int NYS) {
    std::vector<double> f(NY + 1);
    for (int j = 0; j <= NY; ++j) f[j] = std::cos(my * PI * j / NY);
    std::vector<double> Y = dct1_definition(f);
    std::vector<double> Uy(NYS + 1, 0.0);
    for (int k = 0; k <= NYS; ++k) Uy[k] = (k == 0) ? Y[k] / (2.0 * NY) : Y[k] / NY;
    return Uy;
}

// cos(kx·x) 场的期望 x 谱（r2c 半峰：DC=1、内部/Nyquist=0.5，超截断带 0）。
double x_peak(int kxx, int kx, int NX, int NXH, int NXS) {
    if (kx == 0) return (kxx == 0) ? 1.0 : 0.0;
    if (kx > NXS) return 0.0;          // 带外模态不被激发
    return (kxx == kx) ? 0.5 : 0.0;
}
// cos(kz·z) 场的期望 z 谱（c2c 半峰 ±kz，Nyquist 槽亦 0.5，超截断带 0）。
double z_peak(int izz, int kz, int NZ, int NZH, int NZS) {
    if (kz == 0) return (izz == 0) ? 1.0 : 0.0;
    if (std::abs(kz) > NZS) return 0.0;
    return (izz == kz || izz == -kz) ? 0.5 : 0.0;
}

// xyzfft2 差异（§6 表：无 Nyquist 折半/镜像）→ Nyquist 槽峰值 1.0 而非 0.5。
double x2_peak(int kxx, int kx, int NX, int NXH, int NXS) {
    if (kx == 0) return (kxx == 0) ? 1.0 : 0.0;
    if (kx > NXS) return 0.0;
    return (kxx == kx) ? ((kx == NXH) ? 1.0 : 0.5) : 0.0;
}
double z2_peak(int izz, int kz, int NZ, int NZH, int NZS) {
    if (kz == 0) return (izz == 0) ? 1.0 : 0.0;
    if (std::abs(kz) > NZS) return 0.0;
    if (izz != kz && izz != -kz) return 0.0;
    return (std::abs(kz) == NZH) ? 1.0 : 0.5;
}

// 带限随机谱：kx≤NXS、|kz|≤NZS、iy≤NYS。
// 必须是"可实化"谱（实数物理场经 'F' 能产生的谱），否则 'B' 终点 c2r 丢虚部、往返非恒等。
// 往返恒等的充分条件（由单模态探针逐槽验证）：
//   - kx=0 与 kx=NXH（c2r 只取该槽实部）：z 方向 Hermitian —— U(-kz)=conj(U(kz))、kz=0 实、
//     kz=NZH（满谱）实且 ± 相等。理由：c2r 把每条 x 线强制为实，单侧 z 模经 Re() 会折半并镜像。
//   - 其余 kx：z 各槽独立复值（z 的 c2c 链复值精确，kx 系数与 x 共轭槽成对即实）；
//     仅满谱 ±NZH 须相等（xyzfft 求和 / xyzfft2 覆盖都映射到同一 Nyquist 模态）。
std::vector<cdouble> make_band_spectrum(const Mesh& m, uint32_t seed) {
    Spectral s(m);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    const bool fullz = (m.NZS == m.NZH), fullx = (m.NXS == m.NXH);
    for (int iy = 0; iy <= m.NYS; ++iy)
        for (int kx = 0; kx <= m.NXS; ++kx) {
            const bool dc_type = (kx == 0 || (fullx && kx == m.NXH));
            if (dc_type) {
                for (int kz = 0; kz <= m.NZS; ++kz) {
                    if (kz == 0) {
                        s.at(0, iy, kx) = {d(rng), 0.0};
                    } else if (fullz && kz == m.NZH) {
                        const double r = d(rng);
                        s.at(kz, iy, kx) = {r, 0.0};
                        s.at(-kz, iy, kx) = {r, 0.0};
                    } else {
                        const cdouble v{d(rng), d(rng)};
                        s.at(kz, iy, kx) = v;
                        s.at(-kz, iy, kx) = std::conj(v);
                    }
                }
            } else {
                for (int kz = -m.NZS; kz <= m.NZS; ++kz) s.at(kz, iy, kx) = {d(rng), d(rng)};
                if (fullz) s.at(-m.NZH, iy, kx) = s.at(m.NZH, iy, kx);
            }
        }
    return s.v;
}

// ---- naive 2D DFT 参考（xzfft 交叉验证，独立于 FFTW）----
// 按生产语义：x r2c(/NX、满谱 Nyquist×0.5、k>NXS 置零) + z c2c(/NZ、满谱镜像)。
void xzfft_reference(const double* UP, cdouble* Uref, const Mesh& m) {
    const int NX = m.NX, NZ = m.NZ, NXH = m.NXH, NZH = m.NZH;
    const int NXS = m.NXS, NZS = m.NZS, NXSP1 = m.NXSP + 1;
    const bool fullx = (NXS == NXH), fullz = (NZS == NZH);
    const double invNX = 1.0 / NX, invNZ = 1.0 / NZ;
    for (int kxx = 0; kxx <= NXS; ++kxx) {
        for (int kzz = -NZS; kzz <= NZS; ++kzz) {
            if (fullz && kzz == -NZH) continue;  // 镜像后补
            cdouble s{0, 0};
            for (int iz = 0; iz < NZ; ++iz)
                for (int ix = 0; ix < NX; ++ix) {
                    double ph = -2.0 * PI * (static_cast<double>(kxx) * ix / NX +
                                             static_cast<double>(kzz) * iz / NZ);
                    s += UP[static_cast<size_t>(iz) * NX + ix] *
                         cdouble(std::cos(ph), std::sin(ph));
                }
            double norm = invNX * invNZ;
            if (fullx && kxx == NXH) norm *= 0.5;
            if (fullz && kzz == NZH) norm *= 0.5;
            Uref[static_cast<size_t>(kzz + NZS) * NXSP1 + kxx] = s * norm;
        }
        if (fullz)  // 满谱镜像：U(-NZH)=U(+NZH)
            Uref[static_cast<size_t>(-NZH + NZS) * NXSP1 + kxx] =
                Uref[static_cast<size_t>(NZH + NZS) * NXSP1 + kxx];
    }
}

// 报告实际最大偏差（调试用；断言级别见调用处）
double report(const char* what, double md) {
    std::cout << "  [info] " << what << " maxdiff = " << md << '\n';
    return md;
}
}  // namespace

// ===================== F1：channel_xyzfft 'F' 解析（C-FULL / C-TRUNC） =====================
TEST(Xyzfft, ForwardAnalytic) {
    constexpr int kCfgs[2][3] = {{16, 30, 16}, {10, 20, 10}};  // (NXS, NYS, NZS)
    for (int cfg = 0; cfg < 2; ++cfg) {
        const bool full = (cfg == 0);
        const int NXS = kCfgs[cfg][0], NYS = kCfgs[cfg][1], NZS = kCfgs[cfg][2];
        const MeshParm mp = make_meshparm(NXS, NYS, NZS);
        const Mesh mesh(mp, 2.0 * PI, 2.0 * PI, 0, 1);
        Transforms T(mesh, 2.0 * PI, 2.0 * PI, 1.0, "");
        const int NX = mesh.NX, NY = mesh.NY, NZ = mesh.NZ;
        const int NXH = mesh.NXH, NZH = mesh.NZH;

        struct Mode { int kx, my, kz; };
        const std::vector<Mode> modes = full
            ? std::vector<Mode>{{16, 3, 16}, {0, 0, 0}, {5, 7, -9}, {16, 30, 16}}
            : std::vector<Mode>{{3, 3, -4}, {10, 20, 10}, {0, 2, 0}};

        for (const Mode& mo : modes) {
            const int kx = mo.kx, my = mo.my, kz = mo.kz;
            Physical P(mesh);
            for (int iz = 0; iz < NZ; ++iz)
                for (int iy = 0; iy <= mesh.NYP; ++iy)
                    for (int ix = 0; ix < NX; ++ix)
                        P.at(iz, iy, ix) =
                            std::cos(2.0 * PI * kx * ix / NX) *
                            std::cos(my * PI * iy / NY) *
                            std::cos(2.0 * PI * kz * iz / NZ);

            const std::vector<double> Uy = y_spectrum_expected(my, NY, NYS);
            Spectral U(mesh);
            T.channel_xyzfft('F', P.v.data(), U.v.data());

            double md = 0.0;
            for (int iz = -NZS; iz <= NZS; ++iz)
                for (int iy = 0; iy <= NYS; ++iy)
                    for (int kxx = 0; kxx <= NXS; ++kxx) {
                        const cdouble expct = cdouble(x_peak(kxx, kx, NX, NXH, NXS) *
                                                      Uy[iy] * z_peak(iz, kz, NZ, NZH, NZS), 0.0);
                        md = std::max(md, std::abs(U.at(iz, iy, kxx) - expct));
                    }
            EXPECT_LT(report("F1", md), 1e-12)
                << "cfg=" << cfg << " (NXS,NYS,NZS)=(" << NXS << ',' << NYS << ',' << NZS
                << ") mode=(" << kx << ',' << my << ',' << kz << ')';
        }
    }
}

// ===================== F2/F3：channel_xyzfft 'F'↔'B' 往返 =====================
// F3 谱→物理→谱：U0 →'B'→ P →'F'→ U1，断言 U1≈U0。
// F2 物理→谱→物理：P（带限）→'F'→ U →'B'→ P2，断言 P2≈P。
TEST(Xyzfft, Roundtrip) {
    constexpr int kCfgs[2][3] = {{16, 30, 16}, {10, 20, 10}};
    for (int cfg = 0; cfg < 2; ++cfg) {
        const int NXS = kCfgs[cfg][0], NYS = kCfgs[cfg][1], NZS = kCfgs[cfg][2];
        const MeshParm mp = make_meshparm(NXS, NYS, NZS);
        const Mesh mesh(mp, 2.0 * PI, 2.0 * PI, 0, 1);
        Transforms T(mesh, 2.0 * PI, 2.0 * PI, 1.0, "");

        // F3：带限随机谱 → 'B' 合成物理 → 'F' 再提取
        Spectral U0(mesh);
        U0.v = make_band_spectrum(mesh, 0xC0DE + static_cast<uint32_t>(cfg));
        Physical P(mesh);
        T.channel_xyzfft('B', P.v.data(), U0.v.data());
        Spectral U1(mesh);
        T.channel_xyzfft('F', P.v.data(), U1.v.data());
        EXPECT_LT(report("F3", maxdiff(U1, U0)), 1e-12) << "cfg=" << cfg;

        // F2：带限物理场（上面合成的 P）→ 'F' → 'B' 还原
        Spectral U2(mesh);
        T.channel_xyzfft('F', P.v.data(), U2.v.data());
        Physical P2(mesh);
        T.channel_xyzfft('B', P2.v.data(), U2.v.data());
        EXPECT_LT(report("F2", maxdiff(P2, P)), 1e-12) << "cfg=" << cfg;
    }
}

// ===================== F4：xzfft（C-FULL）解析 + naive 2D DFT 交叉 =====================
TEST(Xzfft, AnalyticAndNaiveDft) {
    const int NXS = 16, NYS = 30, NZS = 16;  // C-FULL（满谱 Nyquist 分支）
    const MeshParm mp = make_meshparm(NXS, NYS, NZS);
    const Mesh mesh(mp, 2.0 * PI, 2.0 * PI, 0, 1);
    Transforms T(mesh, 2.0 * PI, 2.0 * PI, 1.0, "");
    const int NX = mesh.NX, NZ = mesh.NZ, NXH = mesh.NXH, NZH = mesh.NZH;

    // 单模态解析
    for (int kz : {3, -5, 0, 16}) {
        for (int kx : {5, 0, 16}) {
            std::vector<double> UP(static_cast<size_t>(NZ) * NX);
            for (int iz = 0; iz < NZ; ++iz)
                for (int ix = 0; ix < NX; ++ix)
                    UP[static_cast<size_t>(iz) * NX + ix] =
                        std::cos(2.0 * PI * kx * ix / NX) * std::cos(2.0 * PI * kz * iz / NZ);
            std::vector<cdouble> U(static_cast<size_t>(mesh.NZS2P1) * (mesh.NXSP + 1), {0, 0});
            T.xzfft(UP.data(), U.data());
            double md = 0.0;
            for (int iz = -NZS; iz <= NZS; ++iz)
                for (int kxx = 0; kxx <= NXS; ++kxx) {
                    const cdouble expct = cdouble(x_peak(kxx, kx, NX, NXH, NXS) *
                                                  z_peak(iz, kz, NZ, NZH, NZS), 0.0);
                    md = std::max(md, std::abs(U[static_cast<size_t>(iz + NZS) * (mesh.NXSP + 1) + kxx] -
                                               expct));
                }
            EXPECT_LT(report("F4-ana", md), 1e-12) << "mode=(" << kx << ',' << kz << ')';
        }
    }

    // 随机场 vs naive 2D DFT
    std::mt19937 rng(0xDA7A);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    std::vector<double> UP(static_cast<size_t>(NZ) * NX);
    for (auto& v : UP) v = d(rng);
    std::vector<cdouble> Uref(static_cast<size_t>(mesh.NZS2P1) * (mesh.NXSP + 1), {0, 0});
    xzfft_reference(UP.data(), Uref.data(), mesh);
    std::vector<cdouble> U(static_cast<size_t>(mesh.NZS2P1) * (mesh.NXSP + 1), {0, 0});
    T.xzfft(UP.data(), U.data());
    double md = 0.0;
    for (size_t i = 0; i < U.size(); ++i) md = std::max(md, std::abs(U[i] - Uref[i]));
    EXPECT_LT(report("F4-dft", md), 1e-12);
}

// ===================== F5：channel_xyzfft2（细网格）解析 + 往返 =====================
// xyzfft2 与 xyzfft 差异（§6 表）：x/z 无 Nyquist 折半/镜像；y 用 NY2。
// 本配置 NX2=NY2=NZ2=NX=NY=NZ，故期望因子：x/z cos 半峰同值，y 归一化用 NY2。
TEST(Xyzfft2, ForwardAndRoundtrip) {
    const int NXS = 16, NYS = 30, NZS = 16;  // C-FULL（含 Nyquist）
    const MeshParm mp = make_meshparm(NXS, NYS, NZS);
    const Mesh mesh(mp, 2.0 * PI, 2.0 * PI, 0, 1);
    Transforms T(mesh, 2.0 * PI, 2.0 * PI, 1.0, "");
    const int NX = mesh.NX, NZ = mesh.NZ, NY2 = mesh.NY2;
    const int NXH = mesh.NXH, NZH = mesh.NZH;

    // F5a 解析：y 期望用 NY2 归一化
    for (int kz : {3, -9, 16}) {
        for (int kx : {5, 16}) {
            const int my = 3;
            const std::vector<double> Uy = y_spectrum_expected(my, NY2, NYS);
            Physical P(mesh);
            for (int iz = 0; iz < NZ; ++iz)
                for (int iy = 0; iy <= mesh.NYP; ++iy)
                    for (int ix = 0; ix < NX; ++ix)
                        P.at(iz, iy, ix) =
                            std::cos(2.0 * PI * kx * ix / NX) *
                            std::cos(my * PI * iy / NY2) *
                            std::cos(2.0 * PI * kz * iz / NZ);
            Spectral U(mesh);
            T.channel_xyzfft2('F', P.v.data(), U.v.data());
            double md = 0.0;
            for (int iz = -NZS; iz <= NZS; ++iz)
                for (int iy = 0; iy <= NYS; ++iy)
                    for (int kxx = 0; kxx <= NXS; ++kxx) {
                        const cdouble expct = cdouble(x2_peak(kxx, kx, NX, NXH, NXS) *
                                                      Uy[iy] * z2_peak(iz, kz, NZ, NZH, NZS), 0.0);
                        md = std::max(md, std::abs(U.at(iz, iy, kxx) - expct));
                    }
            EXPECT_LT(report("F5-ana", md), 1e-12) << "mode=(" << kx << ',' << kz << ')';
        }
    }

    // F5b 往返：带限随机可实化谱 → 'B' → 'F'
    Spectral U0(mesh);
    U0.v = make_band_spectrum(mesh, 0xF5);
    Physical P(mesh);
    T.channel_xyzfft2('B', P.v.data(), U0.v.data());
    Spectral U1(mesh);
    T.channel_xyzfft2('F', P.v.data(), U1.v.data());
    EXPECT_LT(report("F5-rt", maxdiff(U1, U0)), 1e-12);

    // F5c 物理往返：'F' → 'B'
    Spectral U2(mesh);
    T.channel_xyzfft2('F', P.v.data(), U2.v.data());
    Physical P2(mesh);
    T.channel_xyzfft2('B', P2.v.data(), U2.v.data());
    EXPECT_LT(report("F5-rtp", maxdiff(P2, P)), 1e-12);
}

TEST(Xyzfft2, QueuedTransformsShareWorkspaceAndKeepCapturedParametersAlive) {
    for (int retained : {10, 16}) for (std::size_t offset : {0u, 2u}) {
        const auto mp = make_meshparm(retained, retained == 16 ? 30 : 20, retained);
        const Mesh mesh(mp, 2.0 * PI, 2.0 * PI, 0, 1);
        Transforms transforms(mesh, 2.0 * PI, 2.0 * PI, 1.0, "");
        const auto source = make_band_spectrum(mesh, 0xAFF1);
        const auto n = static_cast<std::size_t>(mesh.NZ2) * (mesh.NY2P + 1) * mesh.NX2;
        std::vector<double> expected(n), physical(n + 4, -123.0), second(n);
        auto original = source;
        transforms.channel_xyzfft2('B', expected.data(), original.data());
        std::vector<cdouble> expected_spectrum(source.size());
        transforms.channel_xyzfft2('F', expected.data(), expected_spectrum.data());
        const cdouble guard{-7., 3.};
        std::vector<cdouble> input(source.size() * 3, guard), output(source.size() * 3, guard);
        for (std::size_t i = 0; i < source.size(); ++i) input[3 * i + 1] = source[i];
        transforms.enqueue_xyzfft2('B', physical.data() + offset, input.data() + 1, 3);
        transforms.enqueue_xyzfft2('F', physical.data() + offset, output.data() + 2, 3);
        transforms.enqueue_xyzfft2('B', second.data(), output.data() + 2, 3);
        cpu_device().sync();
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(physical[offset + i], expected[i], 1e-12);
            EXPECT_NEAR(second[i], expected[i], 1e-12);
        }
        for (std::size_t i = 0; i < source.size(); ++i) {
            EXPECT_NEAR(std::abs(output[3 * i + 2] - expected_spectrum[i]), 0., 1e-12);
            EXPECT_EQ(output[3 * i], guard);
            EXPECT_EQ(output[3 * i + 1], guard);
            EXPECT_EQ(input[3 * i + 1], source[i]);
        }
        for (std::size_t i = 0; i < offset; ++i) EXPECT_EQ(physical[i], -123.0);
        for (std::size_t i = offset + n; i < physical.size(); ++i) EXPECT_EQ(physical[i], -123.0);
    }
}
