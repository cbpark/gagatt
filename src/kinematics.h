#ifndef SRC_KINEMATICS_H
#define SRC_KINEMATICS_H

#include <cmath>

namespace gagatt {
inline double lambda(double a, double b, double c) {
    const double sqrt_b = std::sqrt(std::max(0.0, b));
    const double sqrt_c = std::sqrt(std::max(0.0, c));
    return (a - std::pow(sqrt_b + sqrt_c, 2)) *
           (a - std::pow(sqrt_b - sqrt_c, 2));
}

inline double lambda12(double a, double b, double c) {
    const double l = lambda(a, b, c);
    return (l <= 0.0) ? 0.0 : std::sqrt(l);
}
}  // namespace gagatt

#endif  // SRC_KINEMATICS_H
