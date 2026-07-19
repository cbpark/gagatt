#ifndef SRC_SPIN_DENSITY_H
#define SRC_SPIN_DENSITY_H

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <Eigen/Core>
#include <cmath>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "amplitude.h"  // for the weighted template
#include "photon.h"

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

    // Construct from kinematics only (uniform 1/4 helicity average).
    SDMatrixCoefficients(double sqrt_s_hat, double cos_th);

    // construct from pre-computed PolarizationCoefficients
    explicit SDMatrixCoefficients(const PolarizationCoefficients &pol);

    SDMatrixCoefficients(double sqrt_s_hat, double cos_th,
                         const LumiWeights &w);

    double c_nn() const { return cc(0, 0); }
    double c_rr() const { return cc(1, 1); }
    double c_kk() const { return cc(2, 2); }
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
inline double getConcurrence(const SDMatrixCoefficients &sdc) {
    const Matrix4cd rho = spinDensityMatrix(sdc);
    return getConcurrence(rho);
}

inline bool isEntangledByConcurrence(const Matrix4cd &rho) {
    // Threshold accounts for numerical noise in the ComplexEigenSolver.
    return getConcurrence(rho) > 1e-12;
}

bool violatesBellInequality(const SDMatrixCoefficients &sdc);

bool isEntangledByD(const SDMatrixCoefficients &sdc);

// Horodecki measure (numeric value, not just bool)
double horodeckiMeasure(const SDMatrixCoefficients &sdc);

// D = (C_nn - |C_rr + C_kk|) / 3
// Entangled when D < -1/3
inline double entanglementMarker(const Eigen::Matrix3d &cij) {
    return (cij(0, 0) - std::abs(cij(1, 1) + cij(2, 2))) / 3.0;
}
inline double entanglementMarker(const SDMatrixCoefficients &sdc) {
    return entanglementMarker(sdc.cc);
}

// Build the 4x4 spin density matrix from B+, B-, and C_ij directly,
// using the same convention as spinDensityMatrix():
//   rho = (1/4)[ I2xI2 + sum_i B+_i (sigma_i x I2)
//                      + sum_j B-_j (I2 x sigma_j)
//                      + sum_{ij} C_ij (sigma_i x sigma_j) ]
Matrix4cd reconstructRho(const Eigen::Vector3d &bp, const Eigen::Vector3d &bm,
                         const Eigen::Matrix3d &cij);
// Convenience overload when B+ = B- = 0 (unpolarized beams or when only
// C_ij is available).
inline Matrix4cd reconstructRho(const Eigen::Matrix3d &cij) {
    return reconstructRho(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(),
                          cij);
}

// Horodecki parameter m12 = two largest eigenvalues of C * C^T.
double m12FromCij(const Eigen::Matrix3d &cij);
}  // namespace gagatt

#endif  // SRC_SPIN_DENSITY_H
