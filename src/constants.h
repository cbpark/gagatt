#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

namespace gagatt {
// fine-structure constant
inline constexpr double ALPHA = 1.0 / 128.0;

// top quark charge in the unit of electric charge
inline constexpr double QTOP = 2.0 / 3.0;

// squared top quark charge in the unit of electric charge
inline constexpr double QTOP2 = QTOP * QTOP;

// top quark mass
inline constexpr double MTOP = 172.5;

// squared top quark mass
inline constexpr double MTOP2 = MTOP * MTOP;

// top quark decay width
inline constexpr double GAMMATOP = 1.4;

// top quark mass * width
inline constexpr double MGAMMATOP = MTOP * GAMMATOP;

// top quark mass * width squared
inline constexpr double MGAMMATOP2 = MGAMMATOP * MGAMMATOP;
}  // namespace gagatt

#endif  // SRC_CONSTANTS_H
