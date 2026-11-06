// FFTW 变换、波数初始化及 MPI 转置。
//
// 串行转置复用同一工作区；MPI 路径在 x 谱块与 y 物理块间显式 pack/unpack。
// plan_many 参数见 step2 §3.2；缓冲布局 [iz][iy][kx]（行主序，kx 最内）。
//
// 对 step2 §3.2 的一处必要修正：y_dct_ 若按原文"howmany=2*(NXSP+1)*NZS2P1 单次执行"，
// 则 iz 平面之间并非 idist=1（平面间距是 2*(NY+1)*(NXSP+1) 个 double），Fourier 变换线会
// 跨平面错位。正确做法：howmany=2*(NXSP+1)（一个 iz 平面内的 re/im 线，起始连续→idist=1），
// 每 iz 平面执行一次（指针偏移一个平面）。NZS2P1 次执行，纯执行次数差异，无数值差异。
#include "fft/transforms.hpp"
#include "core/cpu_device.hpp"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>

#include "core/mpi_context.hpp"
#include "fft/wisdom.hpp"

namespace {
// Drain submitted work before referenced locals unwind if a later submission fails.
class DrainOnFailure {
  public:
    explicit DrainOnFailure(CpuDevice &device) : device_(device), exceptions_(std::uncaught_exceptions()) {}
    ~DrainOnFailure() {
        if (std::uncaught_exceptions() > exceptions_) {
            try { device_.sync(); } catch (...) {} // Preserve the original exception.
        }
    }
  private:
    CpuDevice &device_;
    int exceptions_;
};

#ifndef CHANNEL_DNS_PURE_MPI
void initialize_fftw_device() {
    struct Runtime {
        Runtime() {
            if (!fftw_init_threads())
                throw std::runtime_error("FFTW thread initialization failed");
            fftw_threads_set_callback(
                [](void *(*work)(char *), char *data, size_t stride, int jobs, void *) {
                    cpu_device().for_each(0, jobs, [&](int i) { work(data + stride * i); });
                },
                nullptr);
        }
        ~Runtime() { fftw_cleanup_threads(); }
    };
    static Runtime runtime;
}
#endif
} // namespace

Transforms::Transforms(const Mesh &mesh, double XL, double ZL, double Re,
                       const std::string &wisdom_file, const MpiContext *mpi)
    : m_(mesh), mpi_(mpi), zoff_(mesh.NZS), LDX_(mesh.NXH + 1), LDX2_(mesh.NX2 / 2 + 1) {
    if (mesh.NPART > 1 && (!mpi_ || mpi_->size() != mesh.NPROC || mpi_->rank() != mesh.my_id))
        throw std::invalid_argument("distributed transforms require the matching MPI context");
#ifndef CHANNEL_DNS_PURE_MPI
    const int threads = cpu_device().threads();
    if (threads > 1)
        initialize_fftw_device();
#endif
    const double PI = std::acos(-1.0);

    // ---- 波数（SETUP 36-49 行）----
    const double ALPHAX = 2.0 * PI / XL;
    const double ALPHAZ = 2.0 * PI / ZL;
    WAVX.resize(mesh.NXSP + 1);
    for (int k = 0; k <= mesh.NXSP; ++k)
        WAVX[k] = (mesh.LSUBX[mesh.INDEX_MPI0] + k) * ALPHAX;
    CWAVX.resize(mesh.NXSP + 1);
    CWAVXRe.resize(mesh.NXSP + 1);
    for (int k = 0; k <= mesh.NXSP; ++k) {
        CWAVX[k] = cdouble(0.0, WAVX[k]);
        CWAVXRe[k] = CWAVX[k] / Re;
    }

    const int nzr = 2 * mesh.NZS + 1;
    WAVZ.assign(nzr, 0.0);
    CWAVZ.assign(nzr, cdouble{0, 0});
    CWAVZRe.assign(nzr, cdouble{0, 0});
    for (int k = -mesh.NZS; k <= mesh.NZS; ++k) {
        const double w = k * ALPHAZ;
        WAVZ[k + zoff_] = w;
        CWAVZ[k + zoff_] = cdouble(0.0, w);
        CWAVZRe[k + zoff_] = cdouble(0.0, w / Re);
    }

    // ---- Chebyshev 权重（SETUP 60-67 行）；0/1 置零（§0.1-4）----
    WL.assign(mesh.NYS + 1, 0.0);
    WD.assign(mesh.NYS + 1, 0.0);
    WR.assign(mesh.NYS + 1, 0.0);
    for (int iy = mesh.NYS; iy >= 2; --iy) {
        WL[iy] = 1.0 / (4.0 * iy * (iy - 1));
        WD[iy] = -1.0 / (2.0 * (iy * iy - 1));
        WR[iy] = 1.0 / (4.0 * iy * (iy + 1));
    }
    WL[2] = 2.0 * WL[2];
    for (int iy = mesh.NYS - 1; iy <= mesh.NYS; ++iy)
        WD[iy] = 0.0;
    for (int iy = mesh.NYS - 3; iy <= mesh.NYS; ++iy)
        WR[iy] = 0.0;

    // ---- 工作缓冲（构造时一次性分配复用，§0.1-7）----
    const size_t NYP1 = mesh.NYP + 1, NY1 = mesh.NY + 1;
    XHAT_.assign(static_cast<size_t>(mesh.NZ) * NYP1 * LDX_, cdouble{0, 0});
    UKZ_.assign(static_cast<size_t>(mesh.NZS2P1) * NYP1 * (mesh.NXS + 1), cdouble{0, 0});
    if (mesh.INDEX_Serial != 1)
        UY_.assign(static_cast<size_t>(mesh.NZS2P1) * NY1 * (mesh.NXSP + 1), cdouble{0, 0});
    UPscratch_.assign(static_cast<size_t>(mesh.NX) * NYP1 * mesh.NZ, 0.0);

    const size_t NY2P1 = mesh.NY2P + 1, NY21 = mesh.NY2 + 1;
    XHAT2_.assign(static_cast<size_t>(mesh.NZ2) * NY2P1 * LDX2_ + 1, cdouble{0, 0});
    UKZ2_.assign(static_cast<size_t>(mesh.NZS2P1) * NY2P1 * (mesh.NXS + 1), cdouble{0, 0});
    if (mesh.INDEX_Serial != 1)
        UY2_.assign(static_cast<size_t>(mesh.NZS2P1) * NY21 * (mesh.NXSP + 1), cdouble{0, 0});
    UP2scratch_.assign(static_cast<size_t>(mesh.NX2) * NY2P1 * mesh.NZ2 + 2, 0.0);

    auto make_transpose_plan = [&](bool fine, bool forward, const char *name) {
        AlltoallvPlan plan;
        plan.send_counts.reserve(mesh.NPART);
        plan.send_displs.reserve(mesh.NPART);
        plan.recv_counts.reserve(mesh.NPART);
        plan.recv_displs.reserve(mesh.NPART);
        const int rank = mesh.INDEX_MPI0;
        const std::size_t local_y = fine ? mesh.y2_count(rank) : mesh.y_count(rank);
        const std::size_t local_x = mesh.x_count(rank);
        std::size_t send_doubles = 0;
        std::size_t recv_doubles = 0;

        auto append = [&](std::vector<int> &counts, std::vector<int> &displs,
                          std::size_t complex_count, std::size_t &total_doubles) {
            if (complex_count > static_cast<std::size_t>(INT_MAX) / 2 ||
                total_doubles > static_cast<std::size_t>(INT_MAX) - 2 * complex_count) {
                throw std::overflow_error(std::string(name) +
                                          " packed buffer exceeds INT_MAX doubles");
            }
            displs.push_back(static_cast<int>(total_doubles));
            counts.push_back(static_cast<int>(2 * complex_count));
            total_doubles += 2 * complex_count;
        };

        for (int peer = 0; peer < mesh.NPART; ++peer) {
            const std::size_t peer_y = fine ? mesh.y2_count(peer) : mesh.y_count(peer);
            const std::size_t peer_x = mesh.x_count(peer);
            const std::size_t forward_send =
                static_cast<std::size_t>(mesh.NZS2P1) * local_y * peer_x;
            const std::size_t forward_recv =
                static_cast<std::size_t>(mesh.NZS2P1) * peer_y * local_x;
            append(plan.send_counts, plan.send_displs, forward ? forward_send : forward_recv,
                   send_doubles);
            append(plan.recv_counts, plan.recv_displs, forward ? forward_recv : forward_send,
                   recv_doubles);
        }
        plan.send_complex = send_doubles / 2;
        plan.recv_complex = recv_doubles / 2;
        return plan;
    };

    coarse_forward_ = make_transpose_plan(false, true, "coarse forward transpose");
    coarse_backward_ = make_transpose_plan(false, false, "coarse backward transpose");
    fine_forward_ = make_transpose_plan(true, true, "fine forward transpose");
    fine_backward_ = make_transpose_plan(true, false, "fine backward transpose");

    assert(coarse_forward_.send_counts == coarse_backward_.recv_counts);
    assert(coarse_forward_.send_displs == coarse_backward_.recv_displs);
    assert(coarse_forward_.recv_counts == coarse_backward_.send_counts);
    assert(coarse_forward_.recv_displs == coarse_backward_.send_displs);
    assert(fine_forward_.send_counts == fine_backward_.recv_counts);
    assert(fine_forward_.send_displs == fine_backward_.recv_displs);
    assert(fine_forward_.recv_counts == fine_backward_.send_counts);
    assert(fine_forward_.recv_displs == fine_backward_.send_displs);
    std::size_t max_send_complex = 0;
    std::size_t max_recv_complex = 0;
    for (const AlltoallvPlan *plan :
         {&coarse_forward_, &coarse_backward_, &fine_forward_, &fine_backward_}) {
        max_send_complex = std::max(max_send_complex, plan->send_complex);
        max_recv_complex = std::max(max_recv_complex, plan->recv_complex);
    }
    send_buffer_.resize(max_send_complex);
    recv_buffer_.resize(max_recv_complex);

    // ---- wisdom 导入（失败不报错，仅打印提示）----
    if (!wisdom_file.empty())
        wisdom::import(wisdom_file);

    // ---- plan 创建（MEASURE 踩缓冲，须在填数据前；缓冲已分配）----
    const int NX = mesh.NX, NZ = mesh.NZ;
    const int NXS1 = mesh.NXS + 1;
    fftw_complex *XH = reinterpret_cast<fftw_complex *>(XHAT_.data());
    fftw_complex *X2 = reinterpret_cast<fftw_complex *>(XHAT2_.data());

#ifndef CHANNEL_DNS_PURE_MPI
    if (threads > 1)
        fftw_plan_with_nthreads(threads);
#endif
    // x：粗。UP[iz][iy][ix] 全体线 → XHAT_[iz][iy][kx]（odist=LDX 为 r2c 必写满 pitch）
    x_r2c_ = fftw_plan_many_dft_r2c(1, &NX, static_cast<int>(NZ * NYP1), UPscratch_.data(), nullptr,
                                    1, NX, XH, nullptr, 1, LDX_, FFTW_MEASURE);
    x_c2r_ = fftw_plan_many_dft_c2r(1, &NX, static_cast<int>(NZ * NYP1), XH, nullptr, 1, LDX_,
                                    UPscratch_.data(), nullptr, 1, NX, FFTW_MEASURE);
    // xzfft 专用：UP[iz][ix] 只有 NZ 条线，输出到 XHAT_ 的 iy=0 板。
    // odist 必须为 NYP1*LDX（iy=0 板的行距），与归一化/z 变换的 iz 步长相一致——
    // 原先误写 LDX 导致 32 行输出挤在开头 32*LDX 内，只有第 0 行落到正确位置，z 谱被摊平。
    x_r2c_xz_ = fftw_plan_many_dft_r2c(1, &NX, NZ, UPscratch_.data(), nullptr, 1, NX, XH, nullptr,
                                       1, static_cast<int>(NYP1 * LDX_), FFTW_MEASURE);

    // z/y plans run concurrently on disjoint planes; no nested FFTW team.
#ifndef CHANNEL_DNS_PURE_MPI
    if (threads > 1)
        fftw_plan_with_nthreads(1);
#endif
    // z：粗，原地。线沿 iz，stride=(NYP+1)*LDX；每 iy 执行一次，howmany=NXS+1（截断谱）
    const int zstride = static_cast<int>(NYP1 * LDX_);
    z_fwd_ = fftw_plan_many_dft(1, &NZ, NXS1, XH, nullptr, zstride, 1, XH, nullptr, zstride, 1,
                                FFTW_FORWARD, FFTW_MEASURE);
    z_bwd_ = fftw_plan_many_dft(1, &NZ, NXS1, XH, nullptr, zstride, 1, XH, nullptr, zstride, 1,
                                FFTW_BACKWARD, FFTW_MEASURE);

    // y：粗，REDFT00 自逆。UY_[iz][iy][kx] 按 double 看：packed slot s=2*kx+parity 为一条线，
    // 线沿 iy（istride=2*(NXSP+1)），同平面内起始连续（idist=1）。每 iz 平面执行一次。
    const int hy = 2 * (mesh.NXSP + 1);
    fftw_r2r_kind kind = FFTW_REDFT00;
    const int NY1i = mesh.NY + 1;
    double *UYd = reinterpret_cast<double *>((mesh.INDEX_Serial == 1 ? UKZ_ : UY_).data());
    y_dct_ = fftw_plan_many_r2r(1, &NY1i, hy, UYd, nullptr, hy, 1, UYd, nullptr, hy, 1, &kind,
                                FFTW_MEASURE);

    // 细网格同构（§6 差异表逐条对应：尺寸族换 NX2/NY2/NZ2，截断谱不变 NXS+1）
    const int NX2 = mesh.NX2, NZ2 = mesh.NZ2;
#ifndef CHANNEL_DNS_PURE_MPI
    if (threads > 1)
        fftw_plan_with_nthreads(threads);
#endif
    x2_r2c_ = fftw_plan_many_dft_r2c(1, &NX2, static_cast<int>(NZ2 * NY2P1), UP2scratch_.data(),
                                     nullptr, 1, NX2, X2, nullptr, 1, LDX2_, FFTW_MEASURE);
    x2_c2r_ = fftw_plan_many_dft_c2r(1, &NX2, static_cast<int>(NZ2 * NY2P1), X2, nullptr, 1, LDX2_,
                                     UP2scratch_.data(), nullptr, 1, NX2, FFTW_MEASURE);
#ifndef CHANNEL_DNS_PURE_MPI
    if (threads > 1)
        fftw_plan_with_nthreads(1);
#endif
    const int z2stride = static_cast<int>(NY2P1 * LDX2_);
    z2_fwd_ = fftw_plan_many_dft(1, &NZ2, NXS1, X2, nullptr, z2stride, 1, X2, nullptr, z2stride, 1,
                                 FFTW_FORWARD, FFTW_MEASURE);
    z2_bwd_ = fftw_plan_many_dft(1, &NZ2, NXS1, X2, nullptr, z2stride, 1, X2, nullptr, z2stride, 1,
                                 FFTW_BACKWARD, FFTW_MEASURE);
    const int NY21i = mesh.NY2 + 1;
    double *UY2d = reinterpret_cast<double *>((mesh.INDEX_Serial == 1 ? UKZ2_ : UY2_).data());
    y2_dct_ = fftw_plan_many_r2r(1, &NY21i, hy, UY2d, nullptr, hy, 1, UY2d, nullptr, hy, 1, &kind,
                                 FFTW_MEASURE);

    if (mesh.INDEX_Serial == 1) {
        // Offset the planning buffers to cover either std::vector alignment class.
        // The extra scalar/complex slots above keep shifted planning in bounds.
        for (int r = 0; r < 2; ++r) {
            double *real = UP2scratch_.data() + 2 * r;
            x2_real_alignment_[r] = fftw_alignment_of(real);
            for (int c = 0; c < 2; ++c) {
                auto *complex = X2 + c;
                x2_complex_alignment_[c] = fftw_alignment_of(reinterpret_cast<double *>(complex));
                const int i = 2 * r + c;
                x2_plane_forward_[i] = fftw_plan_many_dft_r2c(
                    1, &NX2, static_cast<int>(NY2P1), real, nullptr, 1, NX2,
                    complex, nullptr, 1, LDX2_, FFTW_MEASURE);
                x2_plane_backward_[i] = fftw_plan_many_dft_c2r(
                    1, &NX2, static_cast<int>(NY2P1), complex, nullptr, 1, LDX2_,
                    real, nullptr, 1, NX2, FFTW_MEASURE);
            }
        }
    }

    if (!wisdom_file.empty())
        wisdom::export_wisdom(wisdom_file);
}

Transforms::~Transforms() {
    for (auto plan : x2_plane_forward_) if (plan) fftw_destroy_plan(plan);
    for (auto plan : x2_plane_backward_) if (plan) fftw_destroy_plan(plan);
    fftw_destroy_plan(x_r2c_);
    fftw_destroy_plan(x_c2r_);
    fftw_destroy_plan(x_r2c_xz_);
    fftw_destroy_plan(z_fwd_);
    fftw_destroy_plan(z_bwd_);
    fftw_destroy_plan(y_dct_);
    fftw_destroy_plan(x2_r2c_);
    fftw_destroy_plan(x2_c2r_);
    fftw_destroy_plan(z2_fwd_);
    fftw_destroy_plan(z2_bwd_);
    fftw_destroy_plan(y2_dct_);
}

namespace {
static_assert(sizeof(cdouble) == 2 * sizeof(double),
              "MPI transpose requires two-double complex storage");
}

void Transforms::transpose_zx_to_y(const std::vector<cdouble> &UKZ, std::vector<cdouble> &UY,
                                   bool fine) const {
    if (m_.INDEX_Serial == 1) {
        if (m_.NPROC != 1 || UKZ.size() != UY.size())
            throw std::logic_error("invalid serial transpose layout");
        cpu_device().copy_n(UKZ.data(), UKZ.size(), UY.data());
        return;
    }
    const int rank = m_.INDEX_MPI0;
    const int local_y = fine ? m_.y2_count(rank) : m_.y_count(rank);
    const int global_y = fine ? m_.NY2 + 1 : m_.NY + 1;
    const int local_x = m_.x_count(rank);
    const int global_x = m_.NXS + 1;
    const AlltoallvPlan &plan = fine ? fine_forward_ : coarse_forward_;
    cpu_device().fill_n(UY.data(), UY.size(), cdouble{});
    cpu_device().for_each(0, static_cast<std::ptrdiff_t>(m_.NPART) * m_.NZS2P1, [&](int task) {
        const int peer = task / m_.NZS2P1;
        const int iz = task % m_.NZS2P1;
        const int peer_x = m_.x_count(peer);
        std::size_t cursor = static_cast<std::size_t>(plan.send_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz) * (local_y * peer_x);
        for (int iy = 0; iy < local_y; ++iy)
            for (int ix = 0; ix < peer_x; ++ix)
                send_buffer_[cursor++] =
                    UKZ[(static_cast<std::size_t>(iz) * local_y + iy) * global_x +
                        m_.x_begin(peer) + ix];
        assert(cursor == static_cast<std::size_t>(plan.send_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz + 1) * (local_y * peer_x));
    });
    assert(plan.send_complex <= send_buffer_.size());
    assert(plan.recv_complex <= recv_buffer_.size());
    mpi_->alltoallv_doubles(reinterpret_cast<const double *>(send_buffer_.data()), plan.send_counts,
                            plan.send_displs, reinterpret_cast<double *>(recv_buffer_.data()),
                            plan.recv_counts, plan.recv_displs,
                            fine ? "fine forward transpose" : "coarse forward transpose");
    cpu_device().for_each(0, static_cast<std::ptrdiff_t>(m_.NPART) * m_.NZS2P1, [&](int task) {
        const int peer = task / m_.NZS2P1;
        const int iz = task % m_.NZS2P1;
        const int peer_y = fine ? m_.y2_count(peer) : m_.y_count(peer);
        const int peer_y_begin = fine ? m_.y2_begin(peer) : m_.y_begin(peer);
        std::size_t cursor = static_cast<std::size_t>(plan.recv_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz) * (peer_y * local_x);
        for (int iy = 0; iy < peer_y; ++iy)
            for (int ix = 0; ix < local_x; ++ix)
                UY[(static_cast<std::size_t>(iz) * global_y + peer_y_begin + iy) * local_x + ix] =
                    recv_buffer_[cursor++];
        assert(cursor == static_cast<std::size_t>(plan.recv_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz + 1) * (peer_y * local_x));
    });
}
void Transforms::transpose_y_to_zx(const std::vector<cdouble> &UY, std::vector<cdouble> &UKZ,
                                   bool fine) const {
    if (m_.INDEX_Serial == 1) {
        if (m_.NPROC != 1 || UKZ.size() != UY.size())
            throw std::logic_error("invalid serial transpose layout");
        cpu_device().copy_n(UY.data(), UY.size(), UKZ.data());
        return;
    }
    const int rank = m_.INDEX_MPI0;
    const int local_y = fine ? m_.y2_count(rank) : m_.y_count(rank);
    const int global_y = fine ? m_.NY2 + 1 : m_.NY + 1;
    const int local_x = m_.x_count(rank);
    const int global_x = m_.NXS + 1;
    const AlltoallvPlan &plan = fine ? fine_backward_ : coarse_backward_;
    cpu_device().fill_n(UKZ.data(), UKZ.size(), cdouble{});
    cpu_device().for_each(0, static_cast<std::ptrdiff_t>(m_.NPART) * m_.NZS2P1, [&](int task) {
        const int peer = task / m_.NZS2P1;
        const int iz = task % m_.NZS2P1;
        const int peer_y = fine ? m_.y2_count(peer) : m_.y_count(peer);
        const int peer_y_begin = fine ? m_.y2_begin(peer) : m_.y_begin(peer);
        std::size_t cursor = static_cast<std::size_t>(plan.send_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz) * (peer_y * local_x);
        for (int iy = 0; iy < peer_y; ++iy)
            for (int ix = 0; ix < local_x; ++ix)
                send_buffer_[cursor++] =
                    UY[(static_cast<std::size_t>(iz) * global_y + peer_y_begin + iy) * local_x +
                       ix];
        assert(cursor == static_cast<std::size_t>(plan.send_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz + 1) * (peer_y * local_x));
    });
    assert(plan.send_complex <= send_buffer_.size());
    assert(plan.recv_complex <= recv_buffer_.size());
    mpi_->alltoallv_doubles(reinterpret_cast<const double *>(send_buffer_.data()), plan.send_counts,
                            plan.send_displs, reinterpret_cast<double *>(recv_buffer_.data()),
                            plan.recv_counts, plan.recv_displs,
                            fine ? "fine backward transpose" : "coarse backward transpose");
    cpu_device().for_each(0, static_cast<std::ptrdiff_t>(m_.NPART) * m_.NZS2P1, [&](int task) {
        const int peer = task / m_.NZS2P1;
        const int iz = task % m_.NZS2P1;
        const int peer_x = m_.x_count(peer);
        std::size_t cursor = static_cast<std::size_t>(plan.recv_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz) * (local_y * peer_x);
        for (int iy = 0; iy < local_y; ++iy)
            for (int ix = 0; ix < peer_x; ++ix)
                UKZ[(static_cast<std::size_t>(iz) * local_y + iy) * global_x + m_.x_begin(peer) +
                    ix] = recv_buffer_[cursor++];
        assert(cursor == static_cast<std::size_t>(plan.recv_displs[peer]) / 2 +
                             static_cast<std::size_t>(iz + 1) * (local_y * peer_x));
    });
}

// ============================ xzfft（串行，仅正变换） ============================
// UP[iz][ix] → U[iz+zoff][kx]（kx=0..NXSP）。
void Transforms::xzfft(const double *UP, cdouble *U) {
    const int NZ = m_.NZ, NXS = m_.NXS, NZS = m_.NZS, NZH = m_.NZH;
    const int NXH = m_.NXH, NXSP1 = m_.NXSP + 1, NYP1 = m_.NYP + 1;
    const bool fullz = (NZS == NZH), fullx = (NXS == NXH);

    // 1. x：r2c 全 NX 线 → XHAT_[iz][iy=0][kx]；归一化 /NX，Nyquist 满谱 ×0.5，k>NXS 置零
    fftw_execute_dft_r2c(x_r2c_xz_, const_cast<double *>(UP),
                         reinterpret_cast<fftw_complex *>(XHAT_.data()));
    const double invNX = 1.0 / m_.NX;
    cpu_device().for_each(0, NZ, [&](int iz) {
        for (int kx = 0; kx <= NXH; ++kx) {
            cdouble &v = XHAT_[static_cast<size_t>(iz) * NYP1 * LDX_ + kx]; // iy=0 板
            v *= invNX;
            if (kx == NXH && fullx)
                v *= 0.5;
            if (kx > NXS)
                v = 0.0;
        }
    });

    // 2. z：z_fwd_ 在 iy=0 平面上执行一次（howmany=NXS+1），散射进 U（/NZ，满谱镜像）
    fftw_execute_dft(z_fwd_, reinterpret_cast<fftw_complex *>(XHAT_.data()),
                     reinterpret_cast<fftw_complex *>(XHAT_.data()));
    std::fill_n(U, static_cast<size_t>(m_.NZS2P1) * NXSP1, cdouble{0, 0});
    const double invNZ = 1.0 / NZ;
    for (int kx = 0; kx < NXSP1; ++kx) {

        const int global_kx = m_.x_begin(m_.INDEX_MPI0) + kx;
        for (int iz = 0; iz <= NZS; ++iz)
            U[static_cast<size_t>(iz + zoff_) * NXSP1 + kx] =
                XHAT_[static_cast<size_t>(iz) * NYP1 * LDX_ + global_kx] * invNZ;
        for (int iz = -NZS; iz <= -1; ++iz)
            U[static_cast<size_t>(iz + zoff_) * NXSP1 + kx] =
                XHAT_[static_cast<size_t>(NZ + iz) * NYP1 * LDX_ + global_kx] * invNZ;
        if (fullz) {
            U[static_cast<size_t>(NZH + zoff_) * NXSP1 + kx] =
                XHAT_[static_cast<size_t>(NZH) * NYP1 * LDX_ + global_kx] * invNZ * 0.5;
            U[static_cast<size_t>(-NZH + zoff_) * NXSP1 + kx] =
                U[static_cast<size_t>(NZH + zoff_) * NXSP1 + kx]; // 镜像
        }
    }
}

// ==================== channel_xyzfft（粗网格 3D，'F'/'B'） ====================
//
void Transforms::channel_xyzfft(char C, double *UP, cdouble *U, int spectral_stride) {
    if (spectral_stride < 1)
        throw std::invalid_argument("spectral stride must be positive");
    auto &y_work = m_.INDEX_Serial == 1 ? UKZ_ : UY_;
    auto &device = cpu_device();
    const int NX = m_.NX, NY = m_.NY, NZ = m_.NZ;
    const int NYP1 = m_.NYP + 1, NY1 = NY + 1;
    const int NYS = m_.NYS, NXS = m_.NXS, NZS = m_.NZS, NZH = m_.NZH, NXH = m_.NXH;
    const int NXSP1 = m_.NXSP + 1;
    const bool fullz = (NZS == NZH), fullx = (NXS == NXH);
    if (m_.INDEX_Serial == 1 && m_.NPROC != 1)
        throw std::logic_error("invalid serial transpose layout");

    // 布局辅助：XHAT_[iz][iy][kx], UKZ_[iz+zoff][iy][kx], UY_[iz+zoff][iy][kx]
    auto X = [&](int iz, int iy, int kx) -> cdouble & {
        return XHAT_[(static_cast<size_t>(iz) * NYP1 + iy) * LDX_ + kx];
    };
    auto K = [&](int iz, int iy, int kx) -> cdouble & {
        return UKZ_[(static_cast<size_t>(iz + zoff_) * NYP1 + iy) * (NXS + 1) + kx];
    };
    auto Y = [&](int iz, int iy, int kx) -> cdouble & {
        return y_work[(static_cast<size_t>(iz + zoff_) * NY1 + iy) * NXSP1 + kx];
    };
    // U 布局 [iz+zoff][iy][kx]
    auto Ui = [&](int iz, int iy, int kx) -> cdouble & {
        return U[((static_cast<size_t>(iz + zoff_) * (NYS + 1) + iy) * NXSP1 + kx) *
                 spectral_stride];
    };

    if (C == 'F') {
        // 1. x: r2c; normalize retained modes in each z-transform task below.
        fftw_execute_dft_r2c(x_r2c_, UP, reinterpret_cast<fftw_complex *>(XHAT_.data()));
        const double invNX = 1.0 / NX;

        // 2. z：每 iy 执行 z_fwd_，散射进 UKZ_（/NZ，满谱 Nyquist 半+镜像）
        fftw_complex *XH = reinterpret_cast<fftw_complex *>(XHAT_.data());
        const double invNZ = 1.0 / NZ;
        double *UYd = reinterpret_cast<double *>(y_work.data());
        DrainOnFailure drain(device);
        device.enqueue_for_each(0, NYP1, [&](int iy) {
            // Normalize retained x modes immediately before using this y plane.
            for (int iz = 0; iz < NZ; ++iz)
                for (int kx = 0; kx <= NXS; ++kx) {
                    X(iz, iy, kx) *= invNX;
                    if (fullx && kx == NXH)
                        X(iz, iy, kx) *= 0.5;
                }
            fftw_execute_dft(z_fwd_, XH + iy * LDX_, XH + iy * LDX_);
            for (int iz = 0; iz <= NZS; ++iz)
                for (int kx = 0; kx <= NXS; ++kx)
                    K(iz, iy, kx) = X(iz, iy, kx) * invNZ;
            for (int iz = -NZS; iz <= -1; ++iz)
                for (int kx = 0; kx <= NXS; ++kx)
                    K(iz, iy, kx) = X(NZ + iz, iy, kx) * invNZ;
            if (fullz)
                for (int kx = 0; kx <= NXS; ++kx) {
                    K(NZH, iy, kx) = X(NZH, iy, kx) * invNZ * 0.5;
                    K(-NZH, iy, kx) = K(NZH, iy, kx);
                }
        });
        // FIFO batches order the aliased serial workspace without waking the owner.
        if (m_.INDEX_Serial != 1) device.sync();

        // 3. 转置：串行复用工作区，分布式交换
        transpose_zx_to_y(UKZ_, y_work, false);

        // 4. y：每 iz 平面 y_dct_，归一化 UY_[0]/(2NY)、[1..NYS]/NY
        device.enqueue_for_each(-NZS, NZS + 1, [&](int iz) {
            fftw_execute_r2r(y_dct_, UYd + 2 * static_cast<size_t>(iz + zoff_) * NY1 * NXSP1,
                             UYd + 2 * static_cast<size_t>(iz + zoff_) * NY1 * NXSP1);

            for (int iy = 0; iy <= NYS; ++iy) {
                const double scale = iy == 0 ? 1.0 / (2.0 * NY) : 1.0 / NY;
                for (int kx = 0; kx < NXSP1; ++kx)
                    Ui(iz, iy, kx) = Y(iz, iy, kx) * scale;
            }
        });
        device.sync();
        return;
    }

    // ---- 'B'----
    // 1. U → UY_：iy=0..NYS 直拷，iy>NYS 补零，iy=0 槽位 ×2（UMY(0)+UMY(0)）
    double *UYd = reinterpret_cast<double *>(y_work.data());
    fftw_complex *XH = reinterpret_cast<fftw_complex *>(XHAT_.data());
    DrainOnFailure drain(device);
    device.enqueue_for_each(-NZS, NZS + 1, [&](int iz) {
        for (int iy = NYS + 1; iy < NY1; ++iy)
            for (int kx = 0; kx < NXSP1; ++kx)
                Y(iz, iy, kx) = 0.0;
        for (int iy = 0; iy <= NYS; ++iy)
            for (int kx = 0; kx < NXSP1; ++kx) {

                cdouble v = Ui(iz, iy, kx);
                if (iy == 0)
                    v = v + v;
                Y(iz, iy, kx) = v;
            }

        fftw_execute_r2r(y_dct_, UYd + 2 * static_cast<size_t>(iz + zoff_) * NY1 * NXSP1,
                         UYd + 2 * static_cast<size_t>(iz + zoff_) * NY1 * NXSP1);
    });
    // FIFO batches order the aliased serial workspace without waking the owner.
    if (m_.INDEX_Serial != 1) device.sync();

    // 3. 转置：串行复用工作区，分布式交换
    transpose_y_to_zx(y_work, UKZ_, false);

    // 4. Expand z and scale by 0.5 while loading; clear only unused slots.
    device.enqueue_for_each(0, NYP1, [&](int iy) {
        for (int iz = 0; iz < NZ; ++iz) {
            const bool retained = iz <= NZS || iz >= NZ - NZS;
            if (retained) {
                const int kz = iz <= NZS ? iz : iz - NZ;
                for (int kx = 0; kx <= NXS; ++kx)
                    X(iz, iy, kx) = K(kz, iy, kx) * 0.5;
                if (fullz && iz == NZH)
                    for (int kx = 0; kx <= NXS; ++kx)
                        X(iz, iy, kx) = K(NZH, iy, kx) * 0.5 + K(-NZH, iy, kx) * 0.5;
            }
            for (int kx = retained ? NXS + 1 : 0; kx < LDX_; ++kx)
                X(iz, iy, kx) = 0.0;
        }
        fftw_execute_dft(z_bwd_, XH + iy * LDX_, XH + iy * LDX_);
        if (fullx)
            for (int iz = 0; iz < NZ; ++iz)
                X(iz, iy, NXH) *= 2.0;
    });
    device.sync();
    fftw_execute_dft_c2r(x_c2r_, XH, UP);
}

// ==================== channel_xyzfft2（细网格 3D，'F'/'B'） ====================
// 结构与 xyzfft 同构，差异点严格按 step2 §6 差异表：
//   - x 正向 ELSE：UMX(2:NX+1)/NX2，无折半（索引用粗 NX，非笔误）；逆向无加倍
//   - z 正向 ELSE：iz=0..NZH 与 -NZH..-1 独立取值，无折半无镜像；逆向两槽分别放
void Transforms::channel_xyzfft2(char C, double *UP2, cdouble *U, int spectral_stride) {
    xyzfft2_impl(C, UP2, U, spectral_stride, true);
}

void Transforms::enqueue_xyzfft2(char C, double *UP2, cdouble *U, int spectral_stride) {
    if (m_.INDEX_Serial != 1)
        throw std::logic_error("asynchronous fine FFT requires the serial layout");
    xyzfft2_impl(C, UP2, U, spectral_stride, false);
}

void Transforms::enqueue_x2(bool forward, double *physical) {
    cpu_device().enqueue_for_each(0, m_.NZ2, [this, forward, physical](int z) {
        auto *real = physical + static_cast<std::size_t>(z) * (m_.NY2P + 1) * m_.NX2;
        auto *complex = reinterpret_cast<fftw_complex *>(XHAT2_.data()) +
                        static_cast<std::size_t>(z) * (m_.NY2P + 1) * LDX2_;
        const int ra = fftw_alignment_of(real);
        const int ca = fftw_alignment_of(reinterpret_cast<double *>(complex));
        const int r = ra == x2_real_alignment_[0] ? 0 : 1;
        const int c = ca == x2_complex_alignment_[0] ? 0 : 1;
        if (ra != x2_real_alignment_[r] || ca != x2_complex_alignment_[c])
            throw std::invalid_argument("unsupported fine FFT buffer alignment");
        const int i = 2 * r + c;
        if (forward) fftw_execute_dft_r2c(x2_plane_forward_[i], real, complex);
        else fftw_execute_dft_c2r(x2_plane_backward_[i], complex, real);
    });
}

void Transforms::xyzfft2_impl(char C, double *UP2, cdouble *U, int spectral_stride, bool wait) {
    if (spectral_stride < 1)
        throw std::invalid_argument("spectral stride must be positive");
    auto &y_work = m_.INDEX_Serial == 1 ? UKZ2_ : UY2_;
    auto &device = cpu_device();
    const int NX2 = m_.NX2, NY2 = m_.NY2, NZ2 = m_.NZ2;
    const int NY2P1 = m_.NY2P + 1, NY21 = NY2 + 1;
    const int NYS = m_.NYS, NXS = m_.NXS, NZS = m_.NZS;
    const int NXSP1 = m_.NXSP + 1;
    if (m_.INDEX_Serial == 1 && m_.NPROC != 1)
        throw std::logic_error("invalid serial transpose layout");

    auto X = [this, NY2P1](int iz, int iy, int kx) -> cdouble & {
        return XHAT2_[(static_cast<size_t>(iz) * NY2P1 + iy) * LDX2_ + kx];
    };
    auto K = [this, NY2P1, NXS](int iz, int iy, int kx) -> cdouble & {
        return UKZ2_[(static_cast<size_t>(iz + zoff_) * NY2P1 + iy) * (NXS + 1) + kx];
    };
    auto Y = [data = y_work.data(), offset = zoff_, NY21, NXSP1](int iz, int iy, int kx) -> cdouble & {
        return data[(static_cast<size_t>(iz + offset) * NY21 + iy) * NXSP1 + kx];
    };
    auto Ui = [U, offset = zoff_, NYS, NXSP1, spectral_stride](int iz, int iy, int kx) -> cdouble & {
        return U[((static_cast<size_t>(iz + offset) * (NYS + 1) + iy) * NXSP1 + kx) *
                 spectral_stride];
    };

    if (C == 'F') {
        // 1. x: r2c; normalize retained modes locally before the z transform.
        if (m_.INDEX_Serial == 1) enqueue_x2(true, UP2);
        else fftw_execute_dft_r2c(x2_r2c_, UP2, reinterpret_cast<fftw_complex *>(XHAT2_.data()));
        const double invNX = 1.0 / NX2;

        // 2. z：每 iy 执行 z2_fwd_，散射进 UKZ2_（/NZ2，正负频独立，无折半无镜像）
        fftw_complex *X2 = reinterpret_cast<fftw_complex *>(XHAT2_.data());
        const double invNZ = 1.0 / NZ2;
        double *UY2d = reinterpret_cast<double *>(y_work.data());
        DrainOnFailure drain(device);
        device.enqueue_for_each(0, NY2P1, [=, this](int iy) {
            // Normalize retained x modes immediately before using this y plane.
            for (int iz = 0; iz < NZ2; ++iz)
                for (int kx = 0; kx <= NXS; ++kx) {
                    X(iz, iy, kx) *= invNX;
                }
            fftw_execute_dft(z2_fwd_, X2 + iy * LDX2_, X2 + iy * LDX2_);
            for (int iz = 0; iz <= NZS; ++iz)
                for (int kx = 0; kx <= NXS; ++kx)
                    K(iz, iy, kx) = X(iz, iy, kx) * invNZ;
            for (int iz = -NZS; iz <= -1; ++iz)
                for (int kx = 0; kx <= NXS; ++kx)
                    K(iz, iy, kx) = X(NZ2 + iz, iy, kx) * invNZ;
        });
        // FIFO batches order the aliased serial workspace without waking the owner.
        if (m_.INDEX_Serial != 1) device.sync();

        // 3. 转置（串行复用工作区）
        transpose_zx_to_y(UKZ2_, y_work, true);

        // 4. y：每 iz 平面 y2_dct_，归一化 UY2_[0]/(2·NY2)、[1..NYS]/NY2
        device.enqueue_for_each(-NZS, NZS + 1, [=, this](int iz) {
            fftw_execute_r2r(y2_dct_, UY2d + 2 * static_cast<size_t>(iz + zoff_) * NY21 * NXSP1,
                             UY2d + 2 * static_cast<size_t>(iz + zoff_) * NY21 * NXSP1);

            for (int iy = 0; iy <= NYS; ++iy) {
                const double scale = iy == 0 ? 1.0 / (2.0 * NY2) : 1.0 / NY2;
                for (int kx = 0; kx < NXSP1; ++kx)
                    Ui(iz, iy, kx) = Y(iz, iy, kx) * scale;
            }
        });
        if (wait) device.sync();
        return;
    }

    // ---- 'B'----
    // 1. U → UY2_：iy=0..NYS 直拷，iy>NYS 补零，iy=0 槽位 ×2
    double *UY2d = reinterpret_cast<double *>(y_work.data());
    fftw_complex *X2 = reinterpret_cast<fftw_complex *>(XHAT2_.data());
    DrainOnFailure drain(device);
    device.enqueue_for_each(-NZS, NZS + 1, [=, this](int iz) {
        for (int iy = NYS + 1; iy < NY21; ++iy)
            for (int kx = 0; kx < NXSP1; ++kx)
                Y(iz, iy, kx) = 0.0;
        for (int iy = 0; iy <= NYS; ++iy)
            for (int kx = 0; kx < NXSP1; ++kx) {

                cdouble v = Ui(iz, iy, kx);
                if (iy == 0)
                    v = v + v;
                Y(iz, iy, kx) = v;
            }

        fftw_execute_r2r(y2_dct_, UY2d + 2 * static_cast<size_t>(iz + zoff_) * NY21 * NXSP1,
                         UY2d + 2 * static_cast<size_t>(iz + zoff_) * NY21 * NXSP1);
    });
    // FIFO batches order the aliased serial workspace without waking the owner.
    if (m_.INDEX_Serial != 1) device.sync();

    // 3. 转置（串行复用工作区）
    transpose_y_to_zx(y_work, UKZ2_, true);

    // 4. Expand z and scale by 0.5 while loading; clear only unused slots.
    device.enqueue_for_each(0, NY2P1, [=, this](int iy) {
        for (int iz = 0; iz < NZ2; ++iz) {
            const bool retained = iz <= NZS || iz >= NZ2 - NZS;
            if (retained) {
                const int kz = iz <= NZS ? iz : iz - NZ2;
                for (int kx = 0; kx <= NXS; ++kx)
                    X(iz, iy, kx) = K(kz, iy, kx) * 0.5;
                // At a full z spectrum the negative Nyquist slot wins, as in
                // the original positive-then-negative expansion loops.
                if (iz == NZ2 - NZS && iz <= NZS)
                    for (int kx = 0; kx <= NXS; ++kx)
                        X(iz, iy, kx) = K(-NZS, iy, kx) * 0.5;
            }
            for (int kx = retained ? NXS + 1 : 0; kx < LDX2_; ++kx)
                X(iz, iy, kx) = 0.0;
        }
        fftw_execute_dft(z2_bwd_, X2 + iy * LDX2_, X2 + iy * LDX2_);
    });
    // Serial x planes are another ordered batch; distributed plans stay on the owner.
    if (m_.INDEX_Serial == 1) {
        enqueue_x2(false, UP2);
        if (wait) device.sync();
    } else {
        device.sync();
        fftw_execute_dft_c2r(x2_c2r_, X2, UP2);
    }
}
