// 单方向 y 变换单测（step2 §7.1 Y）：REDFT00（非归一化 DCT-I）。
//   Y1 解析：Lobatto 点上 T_m（f_j=cos(mπj/NY)）按定义式算期望，Y_k 仅在 k=m 峰型
//   Y2 往返：REDFT00∘REDFT00 = 2·NY·恒等（n=NY+1）
#include <gtest/gtest.h>
#include <fftw3.h>

#include <cmath>
#include <random>
#include <vector>


namespace {
constexpr double PI = 3.14159265358979323846;

// REDFT00（DCT-I）定义式（§0.1-5，逐字）：
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
}  // namespace

TEST(FftY, AnalyticLobattoModes) {
    const int NY = 32;                 // 生产配置
    const int n = NY + 1;
    std::vector<double> x(n), out(n);
    fftw_plan p = fftw_plan_r2r_1d(n, x.data(), out.data(), FFTW_REDFT00, FFTW_ESTIMATE);

    for (int m : {0, 3, NY}) {
        for (int j = 0; j < n; ++j) x[j] = std::cos(m * PI * j / NY);
        fftw_execute(p);
        std::vector<double> expct = dct1_definition(x);
        for (int k = 0; k < n; ++k)
            EXPECT_NEAR(out[k], expct[k], 1e-11) << "m=" << m << " k=" << k;
        // 峰型：k≠m 处近零（DCT-I 正交基），m 处取定义式算出的峰值
        for (int k = 0; k < n; ++k) {
            if (k != m) { EXPECT_NEAR(out[k], 0.0, 1e-11) << "m=" << m << " k=" << k; }
        }
    }
    fftw_destroy_plan(p);
}

TEST(FftY, Roundtrip) {
    const int NY = 32;
    const int n = NY + 1;
    std::vector<double> x(n), mid(n), out(n);
    std::mt19937 rng(0x7EED);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (auto& v : x) v = d(rng);

    fftw_plan f = fftw_plan_r2r_1d(n, x.data(), mid.data(), FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_plan b = fftw_plan_r2r_1d(n, mid.data(), out.data(), FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(f);
    fftw_execute(b);
    fftw_destroy_plan(f);
    fftw_destroy_plan(b);

    for (int j = 0; j < n; ++j)
        EXPECT_NEAR(out[j], 2.0 * NY * x[j], 1e-11) << "j=" << j;
}
