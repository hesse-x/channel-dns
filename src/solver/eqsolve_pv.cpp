#include "solver/eqsolve_pv.hpp"
#include "core/cpu_device.hpp"
#include "core/complex_ops.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "fft/transforms.hpp"
#include "solver/vector_ops.hpp"

namespace {
using Complex = std::complex<double>;

template <typename Read, typename Write, typename Ad>
void solve_tau_column(int ny, const Transforms &t_, double g, double a, double b, bool zero_mode,
                      Complex lower, Complex upper, double w, Read read_rhs, Write out, Ad ad) {
    const bool normalize = g == 0.0 && a == 0.0 && zero_mode;
    Complex correction{};
    if (normalize) {
        for (int y = 2; y <= ny; y += 2)
            correction += read_rhs(y) / static_cast<double>(y * y - 1);
        if (b != 0.0)
            correction += (lower + upper) / (2.0 * b);
    }
    const auto rhs_at = [&](int y) { return normalize && y == 0 ? correction : read_rhs(y); };
    out(0) = (lower + upper) / 2.0;
    out(1) = (lower - upper) / 2.0;
    for (int iy = 2; iy <= ny - 2; ++iy)
        out(iy) = t_.WL[iy] * rhs_at(iy - 2) + t_.WD[iy] * rhs_at(iy) + t_.WR[iy] * rhs_at(iy + 2);
    out(ny - 1) = t_.WL[ny - 1] * rhs_at(ny - 3) + t_.WD[ny - 1] * rhs_at(ny - 1);
    out(ny) = t_.WL[ny] * rhs_at(ny - 2) + t_.WD[ny] * rhs_at(ny);
    for (int iy = ny - 2; iy >= 2; --iy)
        out(iy) -= t_.WR[iy] * w * ad(iy + 2) * out(iy + 2);

    double ad01[2]{};
    for (int iy = ny; iy >= ny - 1; --iy) {
        const int parity = iy % 2;
        ad01[parity] = (a == 0.0 && b == 0.0) ? static_cast<double>((1 - parity) * iy * iy + parity)
                                              : a + b * iy * iy;
    }
    for (int iy = ny - 2; iy >= 0; --iy) {
        const int parity = iy % 2;
        const double func = (a == 0.0 && b == 0.0)
                                ? static_cast<double>((1 - parity) * iy * iy + parity)
                                : a + b * iy * iy;
        out(parity) -= ad01[parity] * ad(iy + 2) * out(iy + 2);
        ad01[parity] = func - t_.WL[iy + 2] * w * ad01[parity] * ad(iy + 2);
    }
    if (g == 0.0 && a == 0.0 && zero_mode) {
        out(0) = 0.0;
        ad01[0] = 1.0;
    }
    assert(std::isfinite(ad01[0]) && std::isfinite(ad01[1]) && ad01[0] != 0.0 && ad01[1] != 0.0);
    out(0) /= ad01[0];
    out(1) /= ad01[1];
    for (int iy = 2; iy <= ny; ++iy)
        out(iy) = (out(iy) - t_.WL[iy] * w * out(iy - 2)) * ad(iy);
}

// Keep adjacent x modes together: y recurrences remain ordered, while each
// recurrence step visits contiguous lanes suitable for cache reuse and SIMD.
template <int B, typename Read, typename Write, typename Ad, typename Weight,
          typename Lower, typename Upper>
void solve_tau_tile(int ny, int lanes, const Transforms &t, double g, double a, double b,
                    int zero_lane, Read rhs, Write out, Ad ad, Weight weight,
                    Lower lower, Upper upper) {
    // Give complete tiles a compile-time trip count; only the tail uses lanes.
    auto solve = [&]<bool Full>() {
        const int count = Full ? B : lanes;
        Complex correction[B]{};
        const bool normalize = g == 0.0 && a == 0.0 && zero_lane >= 0;
        if (normalize) {
            for (int y = 2; y <= ny; y += 2)
                correction[zero_lane] += rhs(y, zero_lane) / static_cast<double>(y * y - 1);
            if (b != 0.0)
                correction[zero_lane] += (lower(zero_lane) + upper(zero_lane)) / (2.0 * b);
        }
        auto read = [&](int y, int l) {
            return normalize && y == 0 && l == zero_lane ? correction[l] : rhs(y, l);
        };
        for (int l = 0; l < count; ++l) {
            out(0, l) = (lower(l) + upper(l)) / 2.0;
            out(1, l) = (lower(l) - upper(l)) / 2.0;
        }
        for (int y = 2; y <= ny - 2; ++y)
            for (int l = 0; l < count; ++l)
                out(y, l) = t.WL[y] * read(y - 2, l) + t.WD[y] * read(y, l) + t.WR[y] * read(y + 2, l);
        for (int l = 0; l < count; ++l) {
            out(ny - 1, l) = t.WL[ny - 1] * read(ny - 3, l) + t.WD[ny - 1] * read(ny - 1, l);
            out(ny, l) = t.WL[ny] * read(ny - 2, l) + t.WD[ny] * read(ny, l);
        }
        for (int y = ny - 2; y >= 2; --y)
            for (int l = 0; l < count; ++l)
                out(y, l) -= t.WR[y] * weight(l) * ad(y + 2, l) * out(y + 2, l);
        double ad01[2][B]{};
        for (int y = ny; y >= ny - 1; --y) {
            const int q = y % 2;
            const double v = (a == 0.0 && b == 0.0) ? static_cast<double>((1 - q) * y * y + q)
                                                   : a + b * y * y;
            for (int l = 0; l < count; ++l) ad01[q][l] = v;
        }
        for (int y = ny - 2; y >= 0; --y) {
            const int q = y % 2;
            const double v = (a == 0.0 && b == 0.0) ? static_cast<double>((1 - q) * y * y + q)
                                                   : a + b * y * y;
            for (int l = 0; l < count; ++l) {
                out(q, l) -= ad01[q][l] * ad(y + 2, l) * out(y + 2, l);
                ad01[q][l] = v - t.WL[y + 2] * weight(l) * ad01[q][l] * ad(y + 2, l);
            }
        }
        if (normalize) {
            out(0, zero_lane) = 0.0;
            ad01[0][zero_lane] = 1.0;
        }
        for (int l = 0; l < count; ++l) {
            assert(std::isfinite(ad01[0][l]) && std::isfinite(ad01[1][l]) &&
                   ad01[0][l] != 0.0 && ad01[1][l] != 0.0);
            out(0, l) /= ad01[0][l]; out(1, l) /= ad01[1][l];
        }
        for (int y = 2; y <= ny; ++y)
            for (int l = 0; l < count; ++l)
                out(y, l) = (out(y, l) - t.WL[y] * weight(l) * out(y - 2, l)) * ad(y, l);
    };
    if (lanes == B) solve.template operator()<true>();
    else solve.template operator()<false>();
}

} // namespace

EqSolvePV::EqSolvePV(const Params &params, const Mesh &mesh, const Transforms &transforms)
    : p_(params), m_(mesh), t_(transforms), wn_pressure_(make_real_mode_field(mesh)),
      wn_velocity_(make_real_mode_field(mesh)), ad_pressure_(make_real_scalar_field(mesh)),
      ad_velocity_(make_real_scalar_field(mesh)),
      pressure_influence_{make_scalar_field(mesh), make_scalar_field(mesh)},
      velocity_influence_{make_vector_field(mesh), make_vector_field(mesh)},
      influence_divergence_{make_scalar_boundary(mesh), make_scalar_boundary(mesh)},
      pressure_b_(make_scalar_field(mesh)), velocity_b_(make_vector_field(mesh)),
      sigma0_(make_scalar_boundary(mesh)), gradient_workspace_(make_vector_field(mesh)),
      forcing_workspace_(make_scalar_field(mesh)),
      component_bc_workspace_(make_scalar_boundary(mesh)), u_dpdx_unit_(mesh.NYS + 1),
      w_dpdz_unit_(mesh.NYS + 1) {
    pressure_a_ = p_.parm.a_pre0;
    pressure_bcoef_ = p_.parm.b_pre0;
}

void EqSolvePV::build_wn_ad(RealModeField &wn, RealScalarField &ad, double g) {
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix)
            wn.at_z(iz, 0, ix, 0) = -(t_.WAVX[ix] * t_.WAVX[ix] +
                                      t_.WAVZ[iz + t_.zoff()] * t_.WAVZ[iz + t_.zoff()] + g);
    });
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            const double w = wn.at_z(iz, 0, ix, 0);
            ad.at_z(iz, m_.NYS, ix, 0) = 1.0 / (1.0 + t_.WD[m_.NYS] * w);
            ad.at_z(iz, m_.NYS - 1, ix, 0) = 1.0 / (1.0 + t_.WD[m_.NYS - 1] * w);
            for (int iy = m_.NYS - 2; iy >= 2; --iy) {
                const double fc = t_.WR[iy] * t_.WL[iy + 2];
                const double denominator =
                    1.0 + (t_.WD[iy] - fc * w * ad.at_z(iz, iy + 2, ix, 0)) * w;
                assert(std::isfinite(denominator) && denominator != 0.0);
                ad.at_z(iz, iy, ix, 0) = 1.0 / denominator;
            }
        }
    });
}

void EqSolvePV::solve_tau(ScalarField &out, const ScalarField &rhs, double g, double a, double b,
                          const ScalarBoundary &bc, const RealModeField &wn,
                          const RealScalarField &ad, int component) const {
    ScalarField alias_copy;
    const ScalarField *source = &rhs;
    if (out.data() == rhs.data()) {
        alias_copy = rhs;
        source = &alias_copy;
    }
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int z) {
        for (int x = 0; x <= m_.NXSP; ++x)
            solve_tau_column(
                m_.NYS, t_, g, a, b, m_.owns_zero_x_mode() && z == 0 && x == 0, bc.at_z(z, 0, x, 0),
                bc.at_z(z, 1, x, 0), wn.at_z(z, 0, x, 0),
                [&](int y) { return source->at_z(z, y, x, 0); },
                [&](int y) -> Complex & { return out.at_z(z, y, x, component); },
                [&](int y) { return ad.at_z(z, y, x, 0); });
    });
}

void EqSolvePV::solve_tau2(ScalarField &out, const ScalarField &rhs, double, const Pair &a,
                           const Pair &b, const ScalarBoundary &bc, const RealModeField &wn,
                           const RealScalarField &ad) const {
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            const double w = wn.at_z(iz, 0, ix, 0);
            out.at_z(iz, 0, ix, 0) = bc.at_z(iz, 0, ix, 0);
            out.at_z(iz, 1, ix, 0) = bc.at_z(iz, 1, ix, 0);
            for (int iy = 2; iy <= m_.NYS - 2; ++iy)
                out.at_z(iz, iy, ix, 0) = t_.WL[iy] * rhs.at_z(iz, iy - 2, ix, 0) +
                                          t_.WD[iy] * rhs.at_z(iz, iy, ix, 0) +
                                          t_.WR[iy] * rhs.at_z(iz, iy + 2, ix, 0);
            out.at_z(iz, m_.NYS - 1, ix, 0) = t_.WL[m_.NYS - 1] * rhs.at_z(iz, m_.NYS - 3, ix, 0) +
                                              t_.WD[m_.NYS - 1] * rhs.at_z(iz, m_.NYS - 1, ix, 0);
            out.at_z(iz, m_.NYS, ix, 0) = t_.WL[m_.NYS] * rhs.at_z(iz, m_.NYS - 2, ix, 0) +
                                          t_.WD[m_.NYS] * rhs.at_z(iz, m_.NYS, ix, 0);
            for (int iy = m_.NYS - 2; iy >= 2; --iy)
                out.at_z(iz, iy, ix, 0) -=
                    t_.WR[iy] * w * ad.at_z(iz, iy + 2, ix, 0) * out.at_z(iz, iy + 2, ix, 0);
            double aa[2]{}, bb[2]{};
            for (int iy = m_.NYS; iy >= m_.NYS - 1; --iy) {
                const int q = iy % 2;
                aa[q] = (1 - 2 * q) * (a[q] + b[q] * iy * iy);
                bb[q] = a[1 - q] + b[1 - q] * iy * iy;
            }
            for (int iy = m_.NYS - 2; iy >= 0; --iy) {
                const int q = iy % 2;
                const double fa = (1 - 2 * q) * (a[q] + b[q] * iy * iy);
                const double fb = a[1 - q] + b[1 - q] * iy * iy;
                const Complex v = ad.at_z(iz, iy + 2, ix, 0) * out.at_z(iz, iy + 2, ix, 0);
                out.at_z(iz, q, ix, 0) -= aa[q] * v;
                out.at_z(iz, 1 - q, ix, 0) -= bb[q] * v;
                aa[q] = fa - t_.WL[iy + 2] * w * aa[q] * ad.at_z(iz, iy + 2, ix, 0);
                bb[q] = fb - t_.WL[iy + 2] * w * bb[q] * ad.at_z(iz, iy + 2, ix, 0);
            }
            const double det = aa[0] * aa[1] - bb[0] * bb[1];
            assert(std::isfinite(det) && det != 0.0);
            const Complex c0 =
                (out.at_z(iz, 0, ix, 0) * aa[1] - out.at_z(iz, 1, ix, 0) * bb[1]) / det;
            const Complex c1 =
                (aa[0] * out.at_z(iz, 1, ix, 0) - bb[0] * out.at_z(iz, 0, ix, 0)) / det;
            out.at_z(iz, 0, ix, 0) = c0;
            out.at_z(iz, 1, ix, 0) = c1;
            for (int iy = 2; iy <= m_.NYS; ++iy)
                out.at_z(iz, iy, ix, 0) =
                    (out.at_z(iz, iy, ix, 0) - t_.WL[iy] * w * out.at_z(iz, iy - 2, ix, 0)) *
                    ad.at_z(iz, iy, ix, 0);
        }
    });
}

void EqSolvePV::solve_coupled_base(ScalarField &pressure, const ScalarField &rhs_pressure,
                                   const ScalarBoundary &pressure_bc, VectorField &velocity,
                                   const VectorField &rhs_velocity,
                                   const VectorBoundary &velocity_bc) {
    solve_tau(pressure, rhs_pressure, p_.parm.G_pre, pressure_a_, pressure_bcoef_, pressure_bc,
              wn_pressure_, ad_pressure_);
    auto &grad = gradient_workspace_;
    vector_ops::gradient(grad, pressure, m_, t_);
    for (int c = 0; c < 3; ++c) {
        auto &forcing = forcing_workspace_;
        cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
            for (int wall = 0; wall < 2; ++wall)
                for (int ix = 0; ix <= m_.NXSP; ++ix)
                    component_bc_workspace_.at_z(iz, wall, ix, 0) =
                        velocity_bc.at_z(iz, wall, ix, c);
            for (int iy = 0; iy <= m_.NYS; ++iy)
                for (int ix = 0; ix <= m_.NXSP; ++ix)
                    forcing.at_z(iz, iy, ix, 0) =
                        p_.parm.Re * (grad.at_z(iz, iy, ix, c) - rhs_velocity.at_z(iz, iy, ix, c));
        });
        solve_tau(velocity, forcing, g_v_, p_.parm.a_V[c], p_.parm.b_V[c], component_bc_workspace_,
                  wn_velocity_, ad_velocity_, c);
    }
}

void EqSolvePV::build_influence_coefficients() {
    ScalarField rhs_p = make_scalar_field(m_);
    VectorField rhs_v = make_vector_field(m_);
    VectorBoundary bc_v = make_vector_boundary(m_);
    for (int i = 0; i < 2; ++i) {
        ScalarBoundary bc_p = make_scalar_boundary(m_);
        cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
            for (int ix = 0; ix <= m_.NXSP; ++ix) {
                bc_p.at_z(iz, 0, ix, 0) = static_cast<double>(1 - i);
                bc_p.at_z(iz, 1, ix, 0) = static_cast<double>(i);
            }
        });
        solve_coupled_base(pressure_influence_[i], rhs_p, bc_p, velocity_influence_[i], rhs_v,
                           bc_v);
        vector_ops::divergence_boundary(influence_divergence_[i], velocity_influence_[i], m_, t_);
    }
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            Complex denominator = influence_divergence_[0].at_z(iz, 0, ix, 0) *
                                      influence_divergence_[1].at_z(iz, 1, ix, 0) -
                                  influence_divergence_[0].at_z(iz, 1, ix, 0) *
                                      influence_divergence_[1].at_z(iz, 0, ix, 0);
            bool root = m_.owns_zero_x_mode() && ix == 0 && iz == 0;
            if (root)
                denominator = p_.parm.a_V[1] == 0.0
                                  ? Complex(1.0, 0.0)
                                  : 0.25 / influence_divergence_[0].at_z(0, 0, 0, 0);
            assert(std::isfinite(denominator.real()) && std::isfinite(denominator.imag()) &&
                   denominator != Complex{});
            for (int wall = 0; wall < 2; ++wall)
                for (int i = 0; i < 2; ++i)
                    // The global zero mode retains only its independent wall constraint.
                    influence_divergence_[i].at_z(iz, wall, ix, 0) =
                        root ? (wall == 0 ? denominator : Complex{})
                             : influence_divergence_[i].at_z(iz, wall, ix, 0) / denominator;
        }
    });
}

void EqSolvePV::solve_coupled(ScalarField &pressure, const ScalarField &rhs_pressure,
                              const ScalarBoundary &pressure_bc, VectorField &velocity,
                              const VectorField &rhs_velocity, const VectorBoundary &velocity_bc) {
    solve_modes(pressure, &rhs_pressure, pressure_bc, velocity, rhs_velocity, velocity_bc, false);
}

void EqSolvePV::solve_modes(ScalarField &pressure, const ScalarField *rhs_pressure,
                            const ScalarBoundary &pressure_bc, VectorField &velocity,
                            const VectorField &rhs_velocity, const VectorBoundary &velocity_bc,
                            bool correct_sigma, bool wait) {
    const bool fuse_rhs = rhs_pressure == nullptr;
    // Adjacent x modes share input cache lines; measured against batches 2/4/8.
    // Scratch uses [y][lane] and is private to each worker, without changing field layout.
    static constexpr int B = 8;
    auto &device = cpu_device();
    device.launch(Grid1D{static_cast<std::size_t>(m_.NZS2P1), 1}, [this, &pressure, rhs_pressure, &pressure_bc, &velocity, &rhs_velocity,
                     &velocity_bc, correct_sigma, fuse_rhs](Block1D block) {
        thread_local std::vector<Complex> scratch;
        const int ny = m_.NYS, n = ny + 1, column = n * B;
        scratch.resize(static_cast<std::size_t>(6) * column);
        Complex *pc = scratch.data(), *dy = pc + column, *forcing = dy + column;
        // Once forcing is built, each input component can hold its solved output.
        Complex *rv[3] = {forcing + column, forcing + 2 * column, forcing + 3 * column};
        Complex *vc[3] = {rv[0], rv[1], rv[2]};
        // Pressure RHS is dead before the pressure derivative is constructed.
        Complex *rp = dy;
        for (auto slab = block.begin; slab < block.end; ++slab) {
            const int z = static_cast<int>(slab) - m_.NZS;
            for (int xb = 0; xb <= m_.NXSP; xb += B) {
                const int lanes = std::min(B, m_.NXSP + 1 - xb);
                const int zero_lane = m_.owns_zero_x_mode() && z == 0 && xb == 0 ? 0 : -1;
                for (int y = 0; y <= ny; ++y)
                    for (int l = 0; l < lanes; ++l)
                        for (int c = 0; c < 3; ++c)
                            rv[c][y * B + l] = rhs_velocity.at_z(z, y, xb + l, c);
                if (fuse_rhs) {
                    for (int l = 0; l < lanes; ++l) {
                        rp[ny * B + l] = 0.0;
                        rp[(ny - 1) * B + l] = static_cast<double>(2 * ny) * rv[1][ny * B + l];
                    }
                    for (int y = ny - 2; y >= 0; --y)
                        for (int l = 0; l < lanes; ++l)
                            rp[y * B + l] = rp[(y + 2) * B + l] +
                                static_cast<double>(2 * (y + 1)) * rv[1][(y + 1) * B + l];
                    for (int l = 0; l < lanes; ++l) rp[l] *= 0.5;
                    for (int y = 0; y <= ny; ++y)
                        for (int l = 0; l < lanes; ++l) {
                            const int k = y * B + l;
                            rp[k] =
                                    (Complex{} + fourier_derivative(t_.CWAVX[xb + l].imag(), rv[0][k])) + rp[k] +
                                    fourier_derivative(t_.CWAVZ[z + t_.zoff()].imag(), rv[2][k]);
                        }
                }
                solve_tau_tile<B>(ny, lanes, t_, p_.parm.G_pre, pressure_a_, pressure_bcoef_, zero_lane,
                    [&](int y, int l) { return fuse_rhs ? rp[y * B + l] : rhs_pressure->at_z(z, y, xb + l, 0); },
                    [&](int y, int l) -> Complex & { return pc[y * B + l]; },
                    [&](int y, int l) { return ad_pressure_.at_z(z, y, xb + l, 0); },
                    [&](int l) { return wn_pressure_.at_z(z, 0, xb + l, 0); },
                    [&](int l) { return pressure_bc.at_z(z, 0, xb + l, 0); },
                    [&](int l) { return pressure_bc.at_z(z, 1, xb + l, 0); });
                for (int l = 0; l < lanes; ++l) {
                    dy[ny * B + l] = 0.0;
                    dy[(ny - 1) * B + l] = static_cast<double>(2 * ny) * pc[ny * B + l];
                }
                for (int y = ny - 2; y >= 0; --y)
                    for (int l = 0; l < lanes; ++l)
                        dy[y * B + l] = dy[(y + 2) * B + l] +
                            static_cast<double>(2 * (y + 1)) * pc[(y + 1) * B + l];
                for (int l = 0; l < lanes; ++l) dy[l] *= 0.5;
                Complex rhs_normal_tail[2][B];
                for (int l = 0; l < lanes; ++l) {
                    rhs_normal_tail[0][l] = rv[1][ny * B + l];
                    rhs_normal_tail[1][l] = rv[1][(ny - 1) * B + l];
                }
                for (int c = 0; c < 3; ++c) {
                    for (int y = 0; y <= ny; ++y)
                        for (int l = 0; l < lanes; ++l) {
                            const int k = y * B + l;
                            const Complex gradient =
                                c == 0 ? fourier_derivative(t_.CWAVX[xb + l].imag(), pc[k])
                                : (c == 1 ? dy[k]
                                          : fourier_derivative(t_.CWAVZ[z + t_.zoff()].imag(), pc[k]));
                            forcing[k] = p_.parm.Re * (gradient - rv[c][k]);
                        }
                    solve_tau_tile<B>(ny, lanes, t_, g_v_, p_.parm.a_V[c], p_.parm.b_V[c], zero_lane,
                        [&](int y, int l) { return forcing[y * B + l]; },
                        [&](int y, int l) -> Complex & { return vc[c][y * B + l]; },
                        [&](int y, int l) { return ad_velocity_.at_z(z, y, xb + l, 0); },
                        [&](int l) { return wn_velocity_.at_z(z, 0, xb + l, 0); },
                        [&](int l) { return velocity_bc.at_z(z, 0, xb + l, c); },
                        [&](int l) { return velocity_bc.at_z(z, 1, xb + l, c); });
                }
                Complex div[2][B]{}, c0[B], c1[B], sigma[2][B]{};
                for (int y = 0; y <= ny; ++y)
                    for (int l = 0; l < lanes; ++l) {
                        const int k = y * B + l;
                        const Complex tangential = fourier_derivative(t_.CWAVX[xb + l].imag(), vc[0][k]) +
                            fourier_derivative(t_.CWAVZ[z + t_.zoff()].imag(), vc[2][k]);
                        const Complex normal = static_cast<double>(y * y) * vc[1][k];
                        div[0][l] += tangential + normal;
                        div[1][l] += y % 2 == 0 ? tangential - normal : -tangential + normal;
                    }
                for (int l = 0; l < lanes; ++l) {
                    const int x = xb + l;
                    c0[l] = div[1][l] * influence_divergence_[1].at_z(z, 0, x, 0) -
                            div[0][l] * influence_divergence_[1].at_z(z, 1, x, 0);
                    c1[l] = div[0][l] * influence_divergence_[0].at_z(z, 1, x, 0) -
                            div[1][l] * influence_divergence_[0].at_z(z, 0, x, 0);
                    if (correct_sigma) {
                        // Only these corrected tail values are needed before the final pass.
                        const Complex ptop = pc[ny * B + l] +
                            (c0[l] * pressure_influence_[0].at_z(z, ny, x, 0) +
                             c1[l] * pressure_influence_[1].at_z(z, ny, x, 0));
                        auto vtail = [&](int y) {
                            return vc[1][y * B + l] +
                                (c0[l] * velocity_influence_[0].at_z(z, y, x, 1) +
                                 c1[l] * velocity_influence_[1].at_z(z, y, x, 1));
                        };
                        const double w = wn_velocity_.at_z(z, 0, x, 0);
                        sigma[0][l] = (-w * vtail(ny) - p_.parm.Re * rhs_normal_tail[0][l]) /
                                      sigma0_.at_z(z, 0, x, 0);
                        sigma[1][l] = (-w * vtail(ny - 1) + p_.parm.Re *
                            (2.0 * ny * ptop - rhs_normal_tail[1][l])) / sigma0_.at_z(z, 1, x, 0);
                    }
                }
                // Apply both corrections and write final fields in one local pass.
                for (int y = 0; y <= ny; ++y)
                    for (int l = 0; l < lanes; ++l) {
                        const int x = xb + l, k = y * B + l, q = (y + 1) % 2;
                        Complex p = pc[k] + (c0[l] * pressure_influence_[0].at_z(z, y, x, 0) +
                                            c1[l] * pressure_influence_[1].at_z(z, y, x, 0));
                        if (correct_sigma) p += sigma[q][l] * pressure_b_.at_z(z, y, x, 0);
                        pressure.at_z(z, y, x, 0) = p;
                        for (int c = 0; c < 3; ++c) {
                            Complex v = vc[c][k] + (c0[l] * velocity_influence_[0].at_z(z, y, x, c) +
                                                   c1[l] * velocity_influence_[1].at_z(z, y, x, c));
                            if (correct_sigma)
                                v += sigma[c == 1 ? 1 - q : q][l] * velocity_b_.at_z(z, y, x, c);
                            velocity.at_z(z, y, x, c) = v;
                        }
                    }
            }
        }
    });
    if (wait) device.sync();
}

void EqSolvePV::build_b_problem() {
    ScalarField rhs_p = make_scalar_field(m_);
    rhs_p.fill(Complex{});
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            rhs_p.at_z(iz, 0, ix, 0) = static_cast<double>(m_.NYS - 1);
            for (int iy = 1; iy < m_.NYS; ++iy)
                rhs_p.at_z(iz, iy, ix, 0) = static_cast<double>(2 * (m_.NYS - ((iy + 1) % 2)));
        }
    });
    ScalarBoundary bc_p = make_scalar_boundary(m_);
    VectorField rhs_v = make_vector_field(m_);
    VectorBoundary bc_v = make_vector_boundary(m_);
    solve_coupled(pressure_b_, rhs_p, bc_p, velocity_b_, rhs_v, bc_v);
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            const double w = wn_velocity_.at_z(iz, 0, ix, 0);
            sigma0_.at_z(iz, 0, ix, 0) = p_.parm.Re + w * velocity_b_.at_z(iz, m_.NYS, ix, 1);
            sigma0_.at_z(iz, 1, ix, 0) =
                p_.parm.Re + w * velocity_b_.at_z(iz, m_.NYS - 1, ix, 1) -
                p_.parm.Re * (2.0 * m_.NYS) * pressure_b_.at_z(iz, m_.NYS, ix, 0);
        }
    });
    if (m_.owns_zero_x_mode())
        sigma0_.at_z(0, 1, 0, 0) = 1.0;
}

void EqSolvePV::build_pressure_gradient_responses() {
    if (!m_.owns_zero_x_mode())
        return;
    ScalarField rhs = make_scalar_field(m_), solved = make_scalar_field(m_);
    ScalarBoundary bc = make_scalar_boundary(m_);
    rhs.at_z(0, 0, 0, 0) = p_.parm.Re;
    solve_tau(solved, rhs, g_v_, p_.parm.a_V[0], p_.parm.b_V[0], bc, wn_velocity_, ad_velocity_);
    qx_dpdx_unit_ = vector_ops::flow_rate(solved, m_);
    for (int iy = 0; iy <= m_.NYS; ++iy)
        u_dpdx_unit_[iy] = solved.at_z(0, iy, 0, 0).real();
    solve_tau(solved, rhs, g_v_, p_.parm.a_V[2], p_.parm.b_V[2], bc, wn_velocity_, ad_velocity_);
    qz_dpdz_unit_ = vector_ops::flow_rate(solved, m_);
    for (int iy = 0; iy <= m_.NYS; ++iy)
        w_dpdz_unit_[iy] = solved.at_z(0, iy, 0, 0).real();
}

void EqSolvePV::setup(double g_v) {
    g_v_ = g_v;
    build_wn_ad(wn_pressure_, ad_pressure_, p_.parm.G_pre);
    build_wn_ad(wn_velocity_, ad_velocity_, g_v_);
    // Solver setup temporarily uses value boundary conditions for the pressure
    // influence and B problems, then restores the configured pressure boundary.
    pressure_a_ = 1.0;
    pressure_bcoef_ = 0.0;
    build_influence_coefficients();
    build_b_problem();
    pressure_a_ = p_.parm.a_pre0;
    pressure_bcoef_ = p_.parm.b_pre0;
    if (m_.owns_zero_x_mode())
        build_pressure_gradient_responses();
}

void EqSolvePV::apply_sigma_correction(ScalarField &pressure, VectorField &velocity,
                                       const VectorField &rhs_velocity) {
    cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
        for (int ix = 0; ix <= m_.NXSP; ++ix) {
            const double w = wn_velocity_.at_z(iz, 0, ix, 0);
            Complex sigma[2];
            sigma[0] = (-w * velocity.at_z(iz, m_.NYS, ix, 1) -
                        p_.parm.Re * rhs_velocity.at_z(iz, m_.NYS, ix, 1)) /
                       sigma0_.at_z(iz, 0, ix, 0);
            sigma[1] = (-w * velocity.at_z(iz, m_.NYS - 1, ix, 1) +
                        p_.parm.Re * (2.0 * m_.NYS * pressure.at_z(iz, m_.NYS, ix, 0) -
                                      rhs_velocity.at_z(iz, m_.NYS - 1, ix, 1))) /
                       sigma0_.at_z(iz, 1, ix, 0);
            for (int iy = 0; iy <= m_.NYS; ++iy) {
                const int q = (iy + 1) % 2;
                pressure.at_z(iz, iy, ix, 0) += sigma[q] * pressure_b_.at_z(iz, iy, ix, 0);
                velocity.at_z(iz, iy, ix, 0) += sigma[q] * velocity_b_.at_z(iz, iy, ix, 0);
                velocity.at_z(iz, iy, ix, 1) += sigma[1 - q] * velocity_b_.at_z(iz, iy, ix, 1);
                velocity.at_z(iz, iy, ix, 2) += sigma[q] * velocity_b_.at_z(iz, iy, ix, 2);
            }
        }
    });
}
