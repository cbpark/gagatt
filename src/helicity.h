#ifndef SRC_HELICITY_H
#define SRC_HELICITY_H

namespace gagatt {
enum class Helicity { PLUS = 1, MINUS = -1, UNPOL = 0 };

inline constexpr double toDouble(Helicity pol) {
    switch (pol) {
    case Helicity::PLUS:
        return 1.0;
    case Helicity::MINUS:
        return -1.0;
    default:
        return 0.0;
    }
}

constexpr Helicity operator-(Helicity pol) {
    switch (pol) {
    case Helicity::PLUS:
        return Helicity::MINUS;
    case Helicity::MINUS:
        return Helicity::PLUS;
    default:
        return pol;
    }
}
}  // namespace gagatt

#endif  // SRC_HELICITY_H
