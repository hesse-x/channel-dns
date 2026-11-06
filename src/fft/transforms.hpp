// Fourier–Chebyshev transforms backed by FFTW.
// Serial execution reuses workspaces; distributed execution transposes with MPI_Alltoallv.
#pragma once

#include <complex>
#include <array>
#include <string>
#include <vector>

#include <fftw3.h>

#include "core/mesh.hpp"

using cdouble = std::complex<double>;
class MpiContext;

class Transforms {
  public:
    // XL/ZL 用于波数（SETUP 36/44 行）；Re 用于 CWAVXRE/CWAVZRE（=CWAVX/Re）。
    // wisdom_file 为空串则跳过 wisdom 导入/导出（单测用）。
    // 偏离 step2 §3.1 说明：ctor 增加 XL/ZL —— 它们在 Parm 而非 Mesh 中，WAVX/WAVZ 必需。
    Transforms(const Mesh &mesh, double XL, double ZL, double Re, const std::string &wisdom_file,
               const MpiContext *mpi = nullptr);
    ~Transforms();                           // fftw_destroy_plan 全部
    Transforms(const Transforms &) = delete; // plan 不可拷贝
    Transforms &operator=(const Transforms &) = delete;

    // ---- 波数与权重（SETUP 36-67 行，原样命名）----
    std::vector<double> WAVX;                            // (0:NXSP) 本 rank x 波数
    std::vector<double> WAVZ;                            // (2*NZS+1)，下标 iz+NZS（zoff_）
    std::vector<cdouble> CWAVX, CWAVZ, CWAVXRe, CWAVZRe; // 同上偏移约定
    std::vector<double> WL, WD, WR; // (0:NYS)，0/1 置零

    // ---- 三组变换入口 ----
    // spectral_stride counts complex elements; pass U+component, stride=3 for vector fields.
    // Backward transforms only read U; forward transforms only write the selected component.
    // UP 物理场 [iz][iy][ix]（xzfft 为 [iz][ix]）；U 谱场 [iz+NZS][iy][kx]。
    void xzfft(const double *UP, cdouble *U); // UP[iz][ix] → U(0:NXSP,-NZS:NZS)
    void channel_xyzfft(char C, double *UP, cdouble *U, int spectral_stride = 1); // 'F'/'B'，粗网格
    void channel_xyzfft2(char C, double *UP2, cdouble *U,
                         int spectral_stride = 1); // 'F'/'B'，细网格
    // Serial layout only. Buffers and this object must outlive the next device sync.
    void enqueue_xyzfft2(char C, double *UP2, cdouble *U, int spectral_stride = 1);

    int zoff() const { return zoff_; }

  private:
    struct AlltoallvPlan {
        std::vector<int> send_counts;
        std::vector<int> send_displs;
        std::vector<int> recv_counts;
        std::vector<int> recv_displs;
        std::size_t send_complex = 0;
        std::size_t recv_complex = 0;
    };

    const Mesh &m_;
    const MpiContext *mpi_;
    int zoff_;       // NZS（z 偏移）
    int LDX_, LDX2_; // x r2c 输出 pitch：NXH+1 / NX2/2+1（FFTW 必写满，截断事后清零）

    fftw_plan x_r2c_ = nullptr, x_c2r_ = nullptr;   // x：粗
    fftw_plan x2_r2c_ = nullptr, x2_c2r_ = nullptr; // x：细
    fftw_plan z_fwd_ = nullptr, z_bwd_ = nullptr;   // z：粗
    fftw_plan z2_fwd_ = nullptr, z2_bwd_ = nullptr; // z：细
    fftw_plan y_dct_ = nullptr, y2_dct_ = nullptr;  // y：REDFT00 自逆，正逆共用
    fftw_plan x_r2c_xz_ = nullptr; // xzfft 专用（howmany=NZ，odist=(NYP+1)*LDX 落 iy=0 板）
    // Single-thread plane plans executed by pool workers. Four combinations cover
    // the two 16-byte alignment classes of real/complex std::vector buffers.
    std::array<fftw_plan, 4> x2_plane_forward_{}, x2_plane_backward_{};
    std::array<int, 2> x2_real_alignment_{}, x2_complex_alignment_{};

    // 工作缓冲（成员复用，构造时一次性分配）
    std::vector<cdouble> XHAT_; // 粗：NZ × (NYP+1) × LDX  —— x 输出 / z 原地
    std::vector<cdouble> UKZ_;  // 粗：NZS2P1 × (NYP+1) × (NXS+1) —— z 截断后
    std::vector<cdouble> UY_;   // 粗：NZS2P1 × (NY+1) × (NXSP+1) —— y 全线缓冲
    std::vector<cdouble> XHAT2_, UKZ2_, UY2_; // 细：NZ2×(NY2P+1)×LDX2 等同构
    // plan 创建（FFTW_MEASURE）会踩缓冲，需真实尺寸的工作区；执行时用调用方数组替代
    std::vector<double> UPscratch_, UP2scratch_; // NX×(NYP+1)×NZ / NX2×(NY2P+1)×NZ2
    mutable std::vector<cdouble> send_buffer_, recv_buffer_;
    AlltoallvPlan coarse_forward_, coarse_backward_;
    AlltoallvPlan fine_forward_, fine_backward_;

    void xyzfft2_impl(char C, double *UP2, cdouble *U, int spectral_stride, bool wait);
    void enqueue_x2(bool forward, double *physical);

    void transpose_zx_to_y(const std::vector<cdouble> &UKZ, std::vector<cdouble> &UY,
                           bool fine) const;
    void transpose_y_to_zx(const std::vector<cdouble> &UY, std::vector<cdouble> &UKZ,
                           bool fine) const;
};
