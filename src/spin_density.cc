#include "spin_density.h"
#include "amplitude.h"

namespace gagatt {
SDMatrixCoefficients::SDMatrixCoefficients(double sqrt_s_hat, double cos_th) {
    const auto pol = computePolCoeffs(sqrt_s_hat, cos_th);

    bp_matrix << pol.c6 + pol.c8, -pol.c5 - pol.c7, pol.c2 + pol.c4;
    bm_matrix << pol.c6 - pol.c8, -pol.c5 + pol.c7, -pol.c2 + pol.c4;

    c_matrix << pol.c13 - pol.c15, -pol.c14 - pol.c16, -pol.c10 - pol.c12,
        pol.c14 - pol.c16, pol.c13 + pol.c15, pol.c9 + pol.c11,
        pol.c10 - pol.c12, -pol.c9 + pol.c11, -pol.c1 + pol.c3;

    norm_factor = pol.c1 + pol.c3;
}
}  // namespace gagatt
