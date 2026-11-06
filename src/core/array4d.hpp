// Array4D<T> —— D5：header-only 四维数组容器。
//
// 行主序，第 3 维连续。线性索引为 ((i*n1 + j)*n2 + k)*n3 + l。
// at_z 在第 0 维加 z_offset，将负波数映射到非负存储下标。
#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

#include "core/cpu_device.hpp"

template <typename T>
class Array4D {
public:
    Array4D() = default;
    Array4D(int n0, int n1, int n2, int n3) { resize(n0, n1, n2, n3); }

    void resize(int n0, int n1, int n2, int n3) {
        n_[0] = n0; n_[1] = n1; n_[2] = n2; n_[3] = n3;
        buf_.assign(static_cast<size_t>(n0) * n1 * n2 * n3, T{});
    }

    T& operator()(int i, int j, int k, int l) {
        return buf_[idx(i, j, k, l)];
    }
    const T& operator()(int i, int j, int k, int l) const {
        return buf_[idx(i, j, k, l)];
    }

    // z 偏移访问器（D4）：iz ∈ [-NZS, NZS] → 物理存储 iz + z_offset。
    // 只作用于第 0 维；release 下不做边界检查（性能）。
    void set_z_offset(int off) { z_offset_ = off; }
    int  z_offset() const { return z_offset_; }
    T& at_z(int iz, int j, int k, int l) {
        return buf_[idx(iz + z_offset_, j, k, l)];
    }
    const T& at_z(int iz, int j, int k, int l) const {
        return buf_[idx(iz + z_offset_, j, k, l)];
    }

    void fill(const T& v) { cpu_device().fill_n(buf_.data(), buf_.size(), v); }

    // Reuse an allocated field and copy its data in parallel. Shapes may change.
    void copy_from(const Array4D& source) {
        if (this == &source) return;
        if (buf_.size() != source.buf_.size()) buf_.resize(source.buf_.size());
        cpu_device().copy_n(source.buf_.data(), source.buf_.size(), buf_.data());
        std::copy_n(source.n_, 4, n_);
        z_offset_ = source.z_offset_;
    }

    T*       data() { return buf_.data(); }
    const T* data() const { return buf_.data(); }
    T*       raw() { return buf_.data(); }        // MPI/FFTW/IO 直接拿底层指针（D6/D8）
    const T* raw() const { return buf_.data(); }

    int size() const { return static_cast<int>(buf_.size()); }
    int dim(int d) const { return n_[d]; }

private:
    size_t idx(int i, int j, int k, int l) const {
        assert(i >= 0 && i < n_[0] && j >= 0 && j < n_[1] &&
               k >= 0 && k < n_[2] && l >= 0 && l < n_[3]);
        return (static_cast<size_t>(i) * n_[1] + j) * n_[2] * n_[3] +
               static_cast<size_t>(k) * n_[3] + l;
    }

    std::vector<T> buf_;
    int n_[4] = {0, 0, 0, 0};
    int z_offset_ = 0;
};
