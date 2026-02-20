#ifndef SRC_HELICITY_H
#define SRC_HELICITY_H

namespace gagatt {
enum class Helicity : int { PLUS = 1, MINUS = -1, UNPOL = 0 };

inline constexpr double toDouble(Helicity pol) {
    return static_cast<double>(pol);
}

inline Helicity toHelicity(double val) {
    constexpr double eps = 1.0e-10;
    if (val > eps) {
        return Helicity::PLUS;
    } else if (val < -eps) {
        return Helicity::MINUS;
    } else {
        return Helicity::UNPOL;
    }
}

inline constexpr Helicity operator-(Helicity pol) {
    return static_cast<Helicity>(-static_cast<int>(pol));
}
}  // namespace gagatt

#endif  // SRC_HELICITY_H
