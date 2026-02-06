#include "spin_density.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <vector>
#include "amplitude.h"

#ifdef DEBUG
#include <iostream>
#endif

namespace gagatt {
SDMatrixCoefficients::SDMatrixCoefficients(double sqrt_s_hat, double cos_th) {
    const auto pol = computePolCoeffs(sqrt_s_hat, cos_th);

    norm_factor = pol.c1 + pol.c3;
    const double inv_norm =
        (std::abs(norm_factor) > 1e-15) ? 1.0 / norm_factor : 0.0;

    // Polarization vectors (B+ and B-)
    bp << (pol.c6 + pol.c8), -(pol.c5 + pol.c7), (pol.c2 + pol.c4);
    bm << (pol.c6 - pol.c8), -(pol.c5 - pol.c7), -(pol.c2 - pol.c4);

    bp *= inv_norm;
    bm *= inv_norm;

    // Correlation matrix (C_ij)
    cc << (pol.c13 - pol.c15), -(pol.c14 + pol.c16), -(pol.c10 + pol.c12),
        (pol.c14 - pol.c16), (pol.c13 + pol.c15), (pol.c9 + pol.c11),
        (pol.c10 - pol.c12), -(pol.c9 - pol.c11), -(pol.c1 - pol.c3);

    cc *= inv_norm;
}

Matrix4cd spinDensityMatrix(const SDMatrixCoefficients &sdc) {
    using namespace Basis;

    // rho = 1/4 ( I + Bp.sigma*I + I*Bm.sigma + C_ij sigma_i*sigma_j )
    Matrix4cd rho = I2I2;

    // Vector polarizations
    rho.noalias() += sdc.bp(0) * S1I2 + sdc.bp(1) * S2I2 + sdc.bp(2) * S3I2;
    rho.noalias() += sdc.bm(0) * I2S1 + sdc.bm(1) * I2S2 + sdc.bm(2) * I2S3;

    // Spin correlations
    rho.noalias() +=
        sdc.cc(0, 0) * S1S1 + sdc.cc(0, 1) * S1S2 + sdc.cc(0, 2) * S1S3;
    rho.noalias() +=
        sdc.cc(1, 0) * S2S1 + sdc.cc(1, 1) * S2S2 + sdc.cc(1, 2) * S2S3;
    rho.noalias() +=
        sdc.cc(2, 0) * S3S1 + sdc.cc(2, 1) * S3S2 + sdc.cc(2, 2) * S3S3;

    return 0.25 * rho;
}

Matrix4cd partialTransposeB(const Matrix4cd &rho) {
    Matrix4cd pt = rho;

    // Transpose each 2x2 sub-block
    pt.block<2, 2>(0, 0).transposeInPlace();
    pt.block<2, 2>(0, 2).transposeInPlace();
    pt.block<2, 2>(2, 0).transposeInPlace();
    pt.block<2, 2>(2, 2).transposeInPlace();

    return pt;
}

bool isEntangled_PH(const Matrix4cd &rho) {
    Matrix4cd rhoPT = partialTransposeB(rho);

    Eigen::SelfAdjointEigenSolver<Matrix4cd> solver(rhoPT);
#ifdef DEBUG
    std::cerr << "isEntangled_PH:\n " << solver.eigenvalues() << '\n';
#endif

    return (solver.eigenvalues().array() < -1e-12).any();
}

double getConcurrence(const Eigen::Matrix4cd &rho) {
    Matrix4cd rho_tilde = Basis::S2S2 * rho.conjugate() * Basis::S2S2;
    Matrix4cd R = rho * rho_tilde;

    Eigen::ComplexEigenSolver<Matrix4cd> solver(R);
    auto e_vals = solver.eigenvalues();

    std::vector<double> lambdas;
    lambdas.reserve(4);
    for (int i = 0; i < 4; ++i) {
        // use real part (imaginary part should be numerical noise)
        double val = std::max(0.0, e_vals[i].real());
        lambdas.push_back(std::sqrt(val));
    }

    // sort lambdas in descending order: l1 >= l2 >= l3 >= l4
    std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());

    double c = lambdas[0] - lambdas[1] - lambdas[2] - lambdas[3];
    return std::max(0.0, c);
}
}  // namespace gagatt
