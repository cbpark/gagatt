#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <numbers>

namespace gagatt {
// physical constants
inline constexpr double ALPHA = 1.0 / 128.0;
inline constexpr double QTOP = 2.0 / 3.0;
inline constexpr double QTOP2 = QTOP * QTOP;

// top quark properties
inline constexpr double MTOP = 172.5;
inline constexpr double MTOP2 = MTOP * MTOP;
inline constexpr double TTBARTHRES = 2.0 * MTOP;
inline constexpr double GAMMATOP = 1.4;
inline constexpr double MGAMMATOP = MTOP * GAMMATOP;
inline constexpr double MGAMMATOP2 = MGAMMATOP * MGAMMATOP;

// derived constants
inline constexpr double COUPLING_FACTOR =
    8.0 * std::numbers::pi * ALPHA * QTOP2;
}  // namespace gagatt

#endif  // SRC_CONSTANTS_H
