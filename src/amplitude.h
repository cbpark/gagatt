#ifndef SRC_AMPLITUDE_H
#define SRC_AMPLITUDE_H

#include <complex>

namespace gagatt {
using Amplitude = std::complex<double>;

enum class Polarization { PLUS, MINUS };

inline constexpr double toDouble(Polarization pol) {
    return (pol == Polarization::PLUS) ? 1.0 : -1.0;
}

inline constexpr Polarization operator-(Polarization pol) {
    return (pol == Polarization::PLUS) ? Polarization::MINUS
                                       : Polarization::PLUS;
}

Amplitude onShellAmplitude(const double s_hat, const double cos_th,
                           const Polarization lambda1,
                           const Polarization lambda2,
                           const Polarization sigma1,
                           const Polarization sigma2);
}  // namespace gagatt

#endif  // SRC_AMPLITUDE_H
