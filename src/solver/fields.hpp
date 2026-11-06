#pragma once

#include <array>
#include <cassert>
#include <complex>

#include "core/array4d.hpp"
#include "core/mesh.hpp"

using ScalarField = Array4D<std::complex<double>>;
using VectorField = Array4D<std::complex<double>>;
using ScalarBoundary = Array4D<std::complex<double>>;
using VectorBoundary = Array4D<std::complex<double>>;
using RealScalarField = Array4D<double>;
using RealModeField = Array4D<double>;

template <typename T>
Array4D<T> make_field(const Mesh& m, int ny, int components) {
    Array4D<T> field(m.NZS2P1, ny, m.NXSP + 1, components);
    field.set_z_offset(m.NZS);
    return field;
}

inline ScalarField make_scalar_field(const Mesh& m) {
    return make_field<std::complex<double>>(m, m.NYS + 1, 1);
}
inline VectorField make_vector_field(const Mesh& m) {
    return make_field<std::complex<double>>(m, m.NYS + 1, 3);
}
inline ScalarBoundary make_scalar_boundary(const Mesh& m) {
    return make_field<std::complex<double>>(m, 2, 1);
}
inline VectorBoundary make_vector_boundary(const Mesh& m) {
    return make_field<std::complex<double>>(m, 2, 3);
}
inline RealScalarField make_real_scalar_field(const Mesh& m) {
    return make_field<double>(m, m.NYS + 1, 1);
}
inline RealModeField make_real_mode_field(const Mesh& m) {
    return make_field<double>(m, 1, 1);
}

template <typename T>
void assert_field_shape(const Array4D<T>& f, const Mesh& m, int ny, int components) {
    assert(f.dim(0) == m.NZS2P1 && f.dim(1) == ny &&
           f.dim(2) == m.NXSP + 1 && f.dim(3) == components &&
           f.z_offset() == m.NZS);
}

