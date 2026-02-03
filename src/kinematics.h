#ifndef SRC_KINEMATICS_H
#define SRC_KINEMATICS_H

#include <cmath>

namespace gagatt {
inline constexpr double lambda(double x, double y, double z) {
    return x * x + y * y + z * z - 2.0 * (x * y + y * z + z * x);
}

inline double lambda12(double x, double y, double z) {
    const double l = lambda(x, y, z);
    return (l <= 0.0) ? 0.0 : std::sqrt(l);
}
}  // namespace gagatt

#endif  // SRC_KINEMATICS_H
