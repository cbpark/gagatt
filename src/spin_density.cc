#include "spin_density.h"
#include <Eigen/Dense>
#include "amplitude.h"

namespace gagatt {
SpinDensityMatrix spinDensityMatrix(double sqrt_s_hat, double cos_th) {
    const auto pol = getPolCoeffs(sqrt_s_hat, cos_th);

    Eigen::Vector3d bp(pol.c6 + pol.c8, -pol.c5 - pol.c7, pol.c2 + pol.c4);
    Eigen::Vector3d bm(pol.c6 - pol.c8, -pol.c5 + pol.c7,
                            -pol.c2 + pol.c4);
    Eigen::Matrix3d c;
    c << pol.c13 - pol.c15, -pol.c14 - pol.c16, -pol.c10 - pol.c12,
        pol.c14 - pol.c16, pol.c13 + pol.c15, pol.c9 + pol.c11,
        pol.c10 - pol.c12, -pol.c9 + pol.c11, -pol.c1 + pol.c3;

    double norm_factor = pol.c1 + pol.c3;

    return {bp, bm, c, norm_factor};
}
}  // namespace gagatt
