#include "spin_density.h"
#include "amplitude.h"

namespace gagatt {
SDMatrixCoefficients::SDMatrixCoefficients(double sqrt_s_hat, double cos_th) {
    const auto pol = computePolCoeffs(sqrt_s_hat, cos_th);

    norm_factor = pol.c1 + pol.c3;

    bp << pol.c6 + pol.c8, -pol.c5 - pol.c7, pol.c2 + pol.c4;
    bp /= norm_factor;

    bm << pol.c6 - pol.c8, -pol.c5 + pol.c7, -pol.c2 + pol.c4;
    bm /= norm_factor;

    cc << pol.c13 - pol.c15, -pol.c14 - pol.c16, -pol.c10 - pol.c12,
        pol.c14 - pol.c16, pol.c13 + pol.c15, pol.c9 + pol.c11,
        pol.c10 - pol.c12, -pol.c9 + pol.c11, -pol.c1 + pol.c3;
    cc /= norm_factor;
}

Matrix4cd spinDensityMatrix(double sqrt_s_hat, double cos_th) {
    const auto sdc = SDMatrixCoefficients(sqrt_s_hat, cos_th);

    Matrix4cd rho = I2I2;

    rho += sdc.bp(0) * S1I2 + sdc.bp(1) * S2I2 + sdc.bp(2) * S3I2;
    rho += sdc.bm(0) * I2S1 + sdc.bm(1) * I2S2 + sdc.bm(2) * I2S3;

    rho += sdc.cc(0, 0) * S1S1 + sdc.cc(0, 1) * S1S2 + sdc.cc(0, 2) * S1S3;
    rho += sdc.cc(1, 0) * S2S1 + sdc.cc(1, 1) * S2S2 + sdc.cc(1, 2) * S2S3;
    rho += sdc.cc(2, 0) * S3S1 + sdc.cc(2, 1) * S3S2 + sdc.cc(2, 2) * S3S3;

    return 0.25 * rho;
}
}  // namespace gagatt
