#ifndef SRC_SPIN_DENSITY_H
#define SRC_SPIN_DENSITY_H

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <Eigen/Core>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "amplitude.h"  // for the weighted template

namespace gagatt {
using Matrix2cd = Eigen::Matrix2cd;
using Matrix4cd = Eigen::Matrix<std::complex<double>, 4, 4>;

namespace Basis {
inline const Matrix2cd I2 = Matrix2cd::Identity();
inline const Matrix2cd S1 = (Matrix2cd() << 0, 1, 1, 0).finished();
inline const Matrix2cd S2 = (Matrix2cd() << 0, std::complex<double>(0, -1),
                             std::complex<double>(0, 1), 0)
                                .finished();
inline const Matrix2cd S3 = (Matrix2cd() << 1, 0, 0, -1).finished();

inline const Matrix4cd I2I2 = Matrix4cd::Identity();

inline const Matrix4cd S1I2 = Eigen::kroneckerProduct(S1, I2).eval();
inline const Matrix4cd S2I2 = Eigen::kroneckerProduct(S2, I2).eval();
inline const Matrix4cd S3I2 = Eigen::kroneckerProduct(S3, I2).eval();

inline const Matrix4cd I2S1 = Eigen::kroneckerProduct(I2, S1).eval();
inline const Matrix4cd I2S2 = Eigen::kroneckerProduct(I2, S2).eval();
inline const Matrix4cd I2S3 = Eigen::kroneckerProduct(I2, S3).eval();

inline const Matrix4cd S1S1 = Eigen::kroneckerProduct(S1, S1).eval();
inline const Matrix4cd S1S2 = Eigen::kroneckerProduct(S1, S2).eval();
inline const Matrix4cd S1S3 = Eigen::kroneckerProduct(S1, S3).eval();
inline const Matrix4cd S2S1 = Eigen::kroneckerProduct(S2, S1).eval();
inline const Matrix4cd S2S2 = Eigen::kroneckerProduct(S2, S2).eval();
inline const Matrix4cd S2S3 = Eigen::kroneckerProduct(S2, S3).eval();
inline const Matrix4cd S3S1 = Eigen::kroneckerProduct(S3, S1).eval();
inline const Matrix4cd S3S2 = Eigen::kroneckerProduct(S3, S2).eval();
inline const Matrix4cd S3S3 = Eigen::kroneckerProduct(S3, S3).eval();
}  // namespace Basis

struct SDMatrixCoefficients {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // The coefficients are normalized by A = C_1 + C_3.
    Eigen::Vector3d bp, bm;
    Eigen::Matrix3d cc;
    double norm_factor;

    // uniform 1/4 average.
    SDMatrixCoefficients(double sqrt_s_hat, double cos_th);

    // construct from pre-computed PolarizationCoefficients
    explicit SDMatrixCoefficients(const PolarizationCoefficients &pol);

    // weighted constructor
    template <typename W>
    SDMatrixCoefficients(double sqrt_s_hat, double cos_th, W &&weight);

    double c_hat_nn() const { return cc(0, 0) * norm_factor; }
    double c_hat_rr() const { return cc(1, 1) * norm_factor; }
    double c_hat_kk() const { return cc(2, 2) * norm_factor; }
};

Matrix4cd spinDensityMatrix(const SDMatrixCoefficients &sdc);

inline Matrix4cd spinDensityMatrix(double sqrt_s_hat, double cos_th) {
    return spinDensityMatrix(SDMatrixCoefficients{sqrt_s_hat, cos_th});
}

template <typename W>
inline Matrix4cd spinDensityMatrix(double sqrt_s_hat, double cos_th,
                                   W &&weight) {
    return spinDensityMatrix(
        SDMatrixCoefficients{sqrt_s_hat, cos_th, std::forward<W>(weight)});
}

Matrix4cd partialTransposeB(const Matrix4cd &rho);

bool isEntangledByPH(const Matrix4cd &rho);

double negativity(const Matrix4cd &rho);

double getConcurrence(const Matrix4cd &rho);

inline bool isEntangledByConcurrence(const Matrix4cd &rho) {
    // eigenvalue solvers for 4x4 complex matrices typically accumulate more
    // numerical noise than 3x3 real symmetric solvers.
    return getConcurrence(rho) > 1e-12;
}

bool violatesBellInequality(const SDMatrixCoefficients &sdc);

bool isEntangledByD(const SDMatrixCoefficients &sdc);

// Horodecki measure (numeric value, not just bool)
double horodeckiMeasure(const SDMatrixCoefficients &sdc);

double entanglementMarker(const SDMatrixCoefficients &sdc);

// weighted SDMatrixCoefficients constructor
template <typename W>
SDMatrixCoefficients::SDMatrixCoefficients(double sqrt_s_hat, double cos_th,
                                           W &&weight)
    : SDMatrixCoefficients(
          computePolCoeffs(sqrt_s_hat, cos_th, std::forward<W>(weight))) {}
}  // namespace gagatt

#endif  // SRC_SPIN_DENSITY_H
