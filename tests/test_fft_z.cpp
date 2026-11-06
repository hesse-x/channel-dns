// 单方向 z 变换单测（step2 §7.1 Z）：
//   Z1 解析：c2c 对 exp(i·2πmj/N) 单峰于 m（含负频 m=-3、Nyquist m=N/2）
//   Z2 往返：BACKWARD(FORWARD(f)) = N·f
#include <gtest/gtest.h>
#include <fftw3.h>

#include <cmath>
#include <complex>
#include <random>
#include <vector>


namespace {
constexpr double PI = 3.14159265358979323846;

}  // namespace

TEST(FftZ, AnalyticModes) {
    const int N = 32;
    std::vector<fftw_complex> in(N), out(N);
    fftw_plan p = fftw_plan_dft_1d(N, in.data(), out.data(), FFTW_FORWARD, FFTW_ESTIMATE);

    for (int m : {3, -3, N / 2}) {
        for (int j = 0; j < N; ++j) {
            double re = std::cos(2.0 * PI * m * j / N);
            double im = std::sin(2.0 * PI * m * j / N);
            in[j][0] = re; in[j][1] = im;
        }
        fftw_execute(p);
        int peak = (m + N) % N;  // 负频 m=-3 → 槽 N-3
        EXPECT_NEAR(out[peak][0], N, 1e-12) << "m=" << m;
        EXPECT_NEAR(out[peak][1], 0.0, 1e-12) << "m=" << m;
        for (int k = 0; k < N; ++k) {
            if (k == peak) continue;
            EXPECT_NEAR(out[k][0], 0.0, 1e-12) << "m=" << m << " k=" << k;
            EXPECT_NEAR(out[k][1], 0.0, 1e-12) << "m=" << m << " k=" << k;
        }
    }
    fftw_destroy_plan(p);
}

TEST(FftZ, Roundtrip) {
    const int N = 32;
    std::vector<fftw_complex> in(N), mid(N), out(N);
    std::mt19937 rng(0x5EED);
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (int j = 0; j < N; ++j) { in[j][0] = d(rng); in[j][1] = d(rng); }

    fftw_plan f = fftw_plan_dft_1d(N, in.data(), mid.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan b = fftw_plan_dft_1d(N, mid.data(), out.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(f);
    fftw_execute(b);
    fftw_destroy_plan(f);
    fftw_destroy_plan(b);

    for (int j = 0; j < N; ++j) {
        EXPECT_NEAR(out[j][0], N * in[j][0], 1e-12) << "j=" << j;
        EXPECT_NEAR(out[j][1], N * in[j][1], 1e-12) << "j=" << j;
    }
}
