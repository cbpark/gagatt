#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <numbers>

namespace gagatt {
// fine-structure constant
inline constexpr double ALPHA = 1.0 / 128.0;

// top quark charge in the unit of electric charge
inline constexpr double QTOP = 2.0 / 3.0;

// squared top quark charge in the unit of electric charge
inline constexpr double QTOP2 = QTOP * QTOP;

// top quark mass
inline constexpr double MTOP = 172.5;

inline constexpr double TTBARTHRES = 2.0 * MTOP;

// squared top quark mass
inline constexpr double MTOP2 = MTOP * MTOP;

// top quark decay width
inline constexpr double GAMMATOP = 1.4;

// top quark mass * width
inline constexpr double MGAMMATOP = MTOP * GAMMATOP;

// top quark mass * width squared
inline constexpr double MGAMMATOP2 = MGAMMATOP * MGAMMATOP;

inline constexpr double COUPLING_FACTOR =
    8.0 * std::numbers::pi * ALPHA * QTOP2;
}  // namespace gagatt

#endif  // SRC_CONSTANTS_H
