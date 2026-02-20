#ifndef SRC_HELICITY_H
#define SRC_HELICITY_H

namespace gagatt {
enum class Helicity : int { PLUS = 1, MINUS = -1, UNPOL = 0 };

inline constexpr double toDouble(Helicity pol) {
    return static_cast<double>(pol);
}

inline constexpr Helicity toHelicity(double val) {
    return val > 0.0 ? Helicity::PLUS : Helicity::MINUS;
}

inline constexpr Helicity operator-(Helicity pol) {
    return static_cast<Helicity>(-static_cast<int>(pol));
}
}  // namespace gagatt

#endif  // SRC_HELICITY_H
