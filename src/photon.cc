#include "photon.h"
#include <cmath>
#include "amplitude.h"

namespace gagatt {
double sigmaC(double x, Helicity pe, Helicity pc) {
    if (x < 1e-12) return 0.0;

    double c0, c1;

    // Taylor expansion for small x to avoid 1/x^2 catastrophic cancellation
    // Analytical limits as x -> 0:
    // sigma_c0 -> (4/3)x - (4/3)x^2 + (26/15)x^3 + ...
    // sigma_c1 -> -(1/3)x^2 + (5/6)x^3 ...
    if (x < 1e-4) {
        c0 = x * (4.0 / 3.0 - x * (4.0 / 3.0 - 26.0 / 15.0 * x));
        c1 = x * x * (-1.0 / 3.0 + 5.0 / 6.0 * x);
    } else {
        const double x2 = x * x;
        const double inv_x = 1.0 / x;
        const double inv_x2 = 1.0 / x2;
        const double log1px = std::log1p(x);
        const double inv_1px = 1.0 / (1.0 + x);
        const double inv_1px2 = inv_1px * inv_1px;

        c0 = (1.0 - 4.0 * inv_x - 8.0 * inv_x2) * log1px +
             (0.5 + 8.0 * inv_x - 0.5 * inv_1px2);

        c1 = (1.0 + 2.0 * inv_x) * log1px - (2.5 - inv_1px + 0.5 * inv_1px2);
    }

    return c0 + toDouble(pe) * toDouble(pc) * c1;
}
}  // namespace gagatt
