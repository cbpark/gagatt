#ifndef SRC_SPIN_DENSITY_H
#define SRC_SPIN_DENSITY_H

#include <Eigen/Core>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std::complex_literals;

namespace gagatt {
using Matrix2cd = Eigen::Matrix2cd;
using Matrix4cd = Eigen::Matrix<std::complex<double>, 4, 4>;

const Matrix2cd I2 = Matrix2cd::Identity();
const Matrix2cd S1 = (Matrix2cd() << 0, 1, 1, 0).finished();
const Matrix2cd S2 = (Matrix2cd() << 0, -1.0i, 1.0i, 0).finished();
const Matrix2cd S3 = (Matrix2cd() << 1, 0, 0, -1).finished();

const Matrix4cd I2I2 = Matrix4cd::Identity();

const Matrix4cd S1I2 = Eigen::kroneckerProduct(S1, I2);
const Matrix4cd S2I2 = Eigen::kroneckerProduct(S2, I2);
const Matrix4cd S3I2 = Eigen::kroneckerProduct(S3, I2);

const Matrix4cd I2S1 = Eigen::kroneckerProduct(I2, S1);
const Matrix4cd I2S2 = Eigen::kroneckerProduct(I2, S2);
const Matrix4cd I2S3 = Eigen::kroneckerProduct(I2, S3);

const Matrix4cd S1S1 = Eigen::kroneckerProduct(S1, S1);
const Matrix4cd S1S2 = Eigen::kroneckerProduct(S1, S2);
const Matrix4cd S1S3 = Eigen::kroneckerProduct(S1, S3);
const Matrix4cd S2S1 = Eigen::kroneckerProduct(S2, S1);
const Matrix4cd S2S2 = Eigen::kroneckerProduct(S2, S2);
const Matrix4cd S2S3 = Eigen::kroneckerProduct(S2, S3);
const Matrix4cd S3S1 = Eigen::kroneckerProduct(S3, S1);
const Matrix4cd S3S2 = Eigen::kroneckerProduct(S3, S2);
const Matrix4cd S3S3 = Eigen::kroneckerProduct(S3, S3);

struct SDMatrixCoefficients {
    Eigen::Vector3d bp_matrix, bm_matrix;
    Eigen::Matrix3d c_matrix;
    double norm_factor;

    SDMatrixCoefficients(double sqrt_s_hat, double cos_th);
};
}  // namespace gagatt

#endif  // SRC_SPIN_DENSITY_H
