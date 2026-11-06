#pragma once

#include <complex>

// Fourier derivative coefficients are purely imaginary: (i*k)*(a+i*b).
inline std::complex<double> fourier_derivative(double wave, std::complex<double> value) {
    return {-wave * value.imag(), wave * value.real()};
}
