// 单方向 x 变换单测（step2 §7.1 X）：
//   X1 解析：cos/sin 单模态 r2c 仅在 m 槽非零，幅值手算
//   X2 往返：c2r(r2c(f)) = N·f
#include <gtest/gtest.h>
#include <fftw3.h>

#include <cmath>
#include <complex>
#include <random>
#include <vector>


namespace {
constexpr double PI = 3.14159265358979323846;

}  // namespace

TEST(FftX, AnalyticModes) {
    const int N = 32;
    std::vector<double> x(N);
    std::vector<fftw_complex> out(N / 2 + 1);
    fftw_plan p = fftw_plan_dft_r2c_1d(N, x.data(), out.data(), FFTW_ESTIMATE);

    for (int m : {1, N / 2 - 1, N / 2}) {
        for (int j = 0; j < N; ++j)
            x[j] = std::cos(2.0 * PI * m * j / N) + 0.3 * std::sin(2.0 * PI * m * j / N);
        fftw_execute(p);
        if (m != N / 2) {
            // out[m] = (N/2)(1 - 0.3i)，其余槽 0（r2c 只存 0..N/2）
            EXPECT_NEAR(out[m][0], N / 2.0, 1e-12);
            EXPECT_NEAR(out[m][1], -0.3 * N / 2.0, 1e-12);
            for (int k = 0; k <= N / 2; ++k) {
                if (k == m) continue;
                EXPECT_NEAR(out[k][0], 0.0, 1e-12);
                EXPECT_NEAR(out[k][1], 0.0, 1e-12);
            }
        } else {
            // Nyquist：f=cos(πj)=(-1)^j，out[N/2] = N（实）
            EXPECT_NEAR(out[N / 2][0], N, 1e-12);
            EXPECT_NEAR(out[N / 2][1], 0.0, 1e-12);
            for (int k = 0; k < N / 2; ++k) {
                EXPECT_NEAR(out[k][0], 0.0, 1e-12);
                EXPECT_NEAR(out[k][1], 0.0, 1e-12);
            }
        }
    }
    fftw_destroy_plan(p);
}

TEST(FftX, Roundtrip) {
    const int N = 32;
    std::vector<double> x(N), y(N);
    std::vector<fftw_complex> mid(N / 2 + 1);
    std::mt19937 rng(0xC0FFEE);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (auto& v : x) v = d(rng);

    fftw_plan f = fftw_plan_dft_r2c_1d(N, x.data(), mid.data(), FFTW_ESTIMATE);
    fftw_plan b = fftw_plan_dft_c2r_1d(N, mid.data(), y.data(), FFTW_ESTIMATE);
    fftw_execute(f);
    fftw_execute(b);
    fftw_destroy_plan(f);
    fftw_destroy_plan(b);

    for (int j = 0; j < N; ++j)
        EXPECT_NEAR(y[j], N * x[j], 1e-12) << "j=" << j;
}
