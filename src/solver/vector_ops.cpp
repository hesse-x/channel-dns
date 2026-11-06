#include "solver/vector_ops.hpp"
#include "core/cpu_device.hpp"
#include "core/complex_ops.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace vector_ops {
namespace {
using Complex = std::complex<double>;

// One complete Chebyshev recurrence stays in one task; columns are independent.
template <typename Read, typename Write> void derivative_y_line(int ny, Read read, Write write) {
    write(ny) = 0.0;
    if (!ny)
        return;
    write(ny - 1) = static_cast<double>(2 * ny) * read(ny);
    for (int y = ny - 2; y >= 0; --y)
        write(y) = write(y + 2) + static_cast<double>(2 * (y + 1)) * read(y + 1);
    write(0) *= 0.5;
}

template <typename Kernel> void spectral_kernel(const Mesh &mesh, Kernel kernel, bool wait = true) {
    auto &device = cpu_device();
    device.launch(Grid1D{static_cast<std::size_t>(mesh.NZS2P1), 1},
                  [kernel, offset = mesh.NZS](Block1D block) {
                      for (auto z = block.begin; z < block.end; ++z)
                          kernel(static_cast<int>(z) - offset);
                  });
    if (wait) device.sync();
}

size_t physical_size(const Mesh &m) { return static_cast<size_t>(m.NZ2) * (m.NY2P + 1) * m.NX2; }
} // namespace

void derivative(int axis, ScalarField &out, const ScalarField &in, const Mesh &m,
                const Transforms &t) {
    assert_field_shape(out, m, m.NYS + 1, 1);
    assert_field_shape(in, m, m.NYS + 1, 1);
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("derivative axis must be 0, 1, or 2");
    if (axis == 0 || axis == 2) {
        cpu_device().for_each(-m.NZS, m.NZS + 1, [&](int iz) {
            for (int iy = 0; iy <= m.NYS; ++iy)
                for (int ix = 0; ix <= m.NXSP; ++ix)
                    out.at_z(iz, iy, ix, 0) =
                        fourier_derivative((axis == 0 ? t.CWAVX[ix] : t.CWAVZ[iz + t.zoff()]).imag(),
                                           in.at_z(iz, iy, ix, 0));
        });
        return;
    }
    out.fill(Complex{});
    if (m.NYS == 0)
        return;
    cpu_device().for_each(-m.NZS, m.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m.NXSP; ++ix) {
            out.at_z(iz, m.NYS, ix, 0) = 0.0;
            out.at_z(iz, m.NYS - 1, ix, 0) =
                static_cast<double>(2 * m.NYS) * in.at_z(iz, m.NYS, ix, 0);
            for (int iy = m.NYS - 2; iy >= 0; --iy)
                out.at_z(iz, iy, ix, 0) =
                    out.at_z(iz, iy + 2, ix, 0) +
                    static_cast<double>(2 * (iy + 1)) * in.at_z(iz, iy + 1, ix, 0);
            out.at_z(iz, 0, ix, 0) *= 0.5;
        }
    });
}

void gradient(VectorField &out, const ScalarField &in, const Mesh &m, const Transforms &t) {
    assert_field_shape(out, m, m.NYS + 1, 3);
    assert_field_shape(in, m, m.NYS + 1, 1);
    spectral_kernel(m, [&](int z) {
        for (int x = 0; x <= m.NXSP; ++x) {
            derivative_y_line(
                m.NYS, [&](int y) { return in.at_z(z, y, x, 0); },
                [&](int y) -> Complex & { return out.at_z(z, y, x, 1); });
            for (int y = 0; y <= m.NYS; ++y) {
                out.at_z(z, y, x, 0) = fourier_derivative(t.CWAVX[x].imag(), in.at_z(z, y, x, 0));
                out.at_z(z, y, x, 2) =
                    fourier_derivative(t.CWAVZ[z + t.zoff()].imag(), in.at_z(z, y, x, 0));
            }
        }
    });
}

void divergence(ScalarField &out, const VectorField &in, const Mesh &m, const Transforms &t) {
    assert_field_shape(out, m, m.NYS + 1, 1);
    assert_field_shape(in, m, m.NYS + 1, 3);
    spectral_kernel(m, [&](int z) {
        for (int x = 0; x <= m.NXSP; ++x) {
            derivative_y_line(
                m.NYS, [&](int y) { return in.at_z(z, y, x, 1); },
                [&](int y) -> Complex & { return out.at_z(z, y, x, 0); });
            for (int y = 0; y <= m.NYS; ++y)
                out.at_z(z, y, x, 0) =
                    (Complex{} + fourier_derivative(t.CWAVX[x].imag(), in.at_z(z, y, x, 0))) +
                    out.at_z(z, y, x, 0) +
                    fourier_derivative(t.CWAVZ[z + t.zoff()].imag(), in.at_z(z, y, x, 2));
        }
    });
}

void curl(VectorField &out, const VectorField &in, const Mesh &m, const Transforms &t) {
    enqueue_curl(out, in, m, t);
    cpu_device().sync();
}

void enqueue_curl(VectorField &out, const VectorField &in, const Mesh &m, const Transforms &t) {
    assert_field_shape(out, m, m.NYS + 1, 3);
    assert_field_shape(in, m, m.NYS + 1, 3);
    spectral_kernel(m, [&](int z) {
        // Keep each y recurrence ordered while visiting adjacent x modes together.
        for (int x = 0; x <= m.NXSP; ++x) {
            out.at_z(z, m.NYS, x, 0) = 0.0;
            out.at_z(z, m.NYS, x, 2) = 0.0;
        }
        if (m.NYS > 0) {
            const double top = static_cast<double>(2 * m.NYS);
            for (int x = 0; x <= m.NXSP; ++x) {
                out.at_z(z, m.NYS - 1, x, 0) = top * in.at_z(z, m.NYS, x, 2);
                out.at_z(z, m.NYS - 1, x, 2) = top * in.at_z(z, m.NYS, x, 0);
            }
            for (int y = m.NYS - 2; y >= 0; --y)
                for (int x = 0; x <= m.NXSP; ++x) {
                    const double scale = static_cast<double>(2 * (y + 1));
                    out.at_z(z, y, x, 0) = out.at_z(z, y + 2, x, 0) +
                                           scale * in.at_z(z, y + 1, x, 2);
                    out.at_z(z, y, x, 2) = out.at_z(z, y + 2, x, 2) +
                                           scale * in.at_z(z, y + 1, x, 0);
                }
            for (int x = 0; x <= m.NXSP; ++x) {
                out.at_z(z, 0, x, 0) *= 0.5;
                out.at_z(z, 0, x, 2) *= 0.5;
            }
        }
        for (int y = 0; y <= m.NYS; ++y)
            for (int x = 0; x <= m.NXSP; ++x) {
                out.at_z(z, y, x, 0) -=
                    fourier_derivative(t.CWAVZ[z + t.zoff()].imag(), in.at_z(z, y, x, 1));
                out.at_z(z, y, x, 1) =
                    fourier_derivative(t.CWAVZ[z + t.zoff()].imag(), in.at_z(z, y, x, 0)) -
                    fourier_derivative(t.CWAVX[x].imag(), in.at_z(z, y, x, 2));
                out.at_z(z, y, x, 2) =
                    fourier_derivative(t.CWAVX[x].imag(), in.at_z(z, y, x, 1)) - out.at_z(z, y, x, 2);
            }
    }, false);
}

void laplacian(ScalarField &out, const ScalarField &in, const Mesh &m, const Transforms &t) {
    VectorField grad = make_vector_field(m);
    gradient(grad, in, m, t);
    divergence(out, grad, m, t);
}

void divergence_boundary(ScalarBoundary &out, const VectorField &in, const Mesh &m,
                         const Transforms &t) {
    spectral_kernel(m, [&](int z) {
        for (int x = 0; x <= m.NXSP; ++x) {
            Complex lower{}, upper{};
            for (int y = 0; y <= m.NYS; ++y) {
                const Complex tangential =
                    fourier_derivative(t.CWAVX[x].imag(), in.at_z(z, y, x, 0)) +
                    fourier_derivative(t.CWAVZ[z + t.zoff()].imag(), in.at_z(z, y, x, 2));
                const Complex normal = static_cast<double>(y * y) * in.at_z(z, y, x, 1);
                lower += tangential + normal;
                upper += (y % 2 == 0) ? tangential - normal : -tangential + normal;
            }
            out.at_z(z, 0, x, 0) = lower;
            out.at_z(z, 1, x, 0) = upper;
        }
    });
}

double max_distance(const ScalarField &lhs, const ScalarField &rhs, const Mesh &) {
    assert(lhs.size() == rhs.size());
    double result = 0.0;
    for (int i = 0; i < lhs.size(); ++i)
        result = std::max(result, std::abs(lhs.data()[i] - rhs.data()[i]));
    return result;
}

double flow_rate(const ScalarField &in, const Mesh &m) {
    double q = 0.0;
    for (int iy = 0; iy <= m.NYS; iy += 2)
        q -= in.at_z(0, iy, 0, 0).real() / static_cast<double>(iy * iy - 1);
    return q;
}

void ProductWorkspace::prepare(const Mesh &mesh) {
    const auto n = physical_size(mesh);
    for (auto &field : a)
        field.resize(n);
    for (auto &field : b)
        field.resize(n);
    product.resize(n);
}

void dot_product(ScalarField &out, const VectorField &a, const VectorField &b, const Mesh &m,
                 Transforms &t, ProductWorkspace *workspace) {
    ProductWorkspace local;
    auto &storage = workspace ? *workspace : local;
    storage.prepare(m);
    auto &ap = storage.a;
    auto &bp = storage.b;
    const size_t n = physical_size(m);
    for (int c = 0; c < 3; ++c) {
        t.channel_xyzfft2('B', ap[c].data(), const_cast<cdouble *>(a.data()) + c, 3);
        t.channel_xyzfft2('B', bp[c].data(), const_cast<cdouble *>(b.data()) + c, 3);
    }
    auto &product = storage.product;
    cpu_device().for_each(
        0, n,
        [&](size_t i) {
            product[i] = ap[0][i] * bp[0][i] + ap[1][i] * bp[1][i] + ap[2][i] * bp[2][i];
        },
        2048);
    t.channel_xyzfft2('F', product.data(), out.data());
}

void cross_product(VectorField &out, const VectorField &a, const VectorField &b, const Mesh &m,
                   Transforms &t, ProductWorkspace *workspace) {
    ProductWorkspace local;
    auto &storage = workspace ? *workspace : local;
    storage.prepare(m);
    auto &ap = storage.a;
    auto &bp = storage.b;
    const size_t n = physical_size(m);
    for (int c = 0; c < 3; ++c) {
        t.channel_xyzfft2('B', ap[c].data(), const_cast<cdouble *>(a.data()) + c, 3);
        t.channel_xyzfft2('B', bp[c].data(), const_cast<cdouble *>(b.data()) + c, 3);
    }
    for (auto &extra : storage.cross_extra) extra.resize(n);
    double *product[3] = {storage.product.data(), storage.cross_extra[0].data(),
                          storage.cross_extra[1].data()};
    // Load all six inputs once, then compute all three outputs in the same pass.
    cpu_device().for_each(0, n, [&](size_t i) {
        const double a0 = ap[0][i], a1 = ap[1][i], a2 = ap[2][i];
        const double b0 = bp[0][i], b1 = bp[1][i], b2 = bp[2][i];
        product[0][i] = a1 * b2 - a2 * b1;
        product[1][i] = a2 * b0 - a0 * b2;
        product[2][i] = a0 * b1 - a1 * b0;
    }, 2048);
    for (int c = 0; c < 3; ++c)
        t.channel_xyzfft2('F', product[c], out.data() + c, 3);
}
void enqueue_cross_product(VectorField &out, const VectorField &a, const VectorField &b,
                           const Mesh &m, Transforms &t, ProductWorkspace &storage) {
    storage.prepare(m);
    const auto n = physical_size(m);
    for (auto &extra : storage.cross_extra) extra.resize(n);
    for (int c = 0; c < 3; ++c) {
        t.enqueue_xyzfft2('B', storage.a[c].data(), const_cast<cdouble *>(a.data()) + c, 3);
        t.enqueue_xyzfft2('B', storage.b[c].data(), const_cast<cdouble *>(b.data()) + c, 3);
    }
    const std::array<const double *, 3> ap{storage.a[0].data(), storage.a[1].data(), storage.a[2].data()};
    const std::array<const double *, 3> bp{storage.b[0].data(), storage.b[1].data(), storage.b[2].data()};
    const std::array<double *, 3> product{storage.product.data(), storage.cross_extra[0].data(),
                                         storage.cross_extra[1].data()};
    cpu_device().enqueue_for_each(0, n, [ap, bp, product](std::size_t i) {
        const double a0 = ap[0][i], a1 = ap[1][i], a2 = ap[2][i];
        const double b0 = bp[0][i], b1 = bp[1][i], b2 = bp[2][i];
        product[0][i] = a1 * b2 - a2 * b1;
        product[1][i] = a2 * b0 - a0 * b2;
        product[2][i] = a0 * b1 - a1 * b0;
    }, 2048);
    for (int c = 0; c < 3; ++c)
        t.enqueue_xyzfft2('F', product[c], out.data() + c, 3);
}

} // namespace vector_ops
