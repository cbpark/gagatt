#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <numbers>

namespace gagatt {
// physical constants
inline constexpr double ALPHA = 1.0 / 128.0;
inline constexpr double QTOP = 2.0 / 3.0;
inline constexpr double QTOP2 = QTOP * QTOP;
inline constexpr int NC = 3;  // color factor

// top quark properties
inline constexpr double MTOP = 172.5;
inline constexpr double MTOP2 = MTOP * MTOP;
inline constexpr double TTBARTHRES = 2.0 * MTOP;
inline constexpr double GAMMATOP = 1.4;
inline constexpr double MGAMMATOP = MTOP * GAMMATOP;
inline constexpr double MGAMMATOP2 = MGAMMATOP * MGAMMATOP;
inline constexpr double BRLL = 0.0455; // BR(t tbar --> l l)
inline constexpr double BRLJ = 0.2877; // BR(t tbar --> l j)

// derived constants
inline constexpr double COUPLING_FACTOR =
    8.0 * std::numbers::pi * ALPHA * QTOP2;

// Conversion: 1 GeV^{-2} = 0.3894 mb = 3.894e11 fb
inline constexpr double GEV2_TO_FB = 3.893793721e11;  // exact to 10 sig figs
}  // namespace gagatt

#endif  // SRC_CONSTANTS_H
