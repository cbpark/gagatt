#ifndef SRC_SPIN_DENSITY_H
#define SRC_SPIN_DENSITY_H

#include <Eigen/Core>

namespace gagatt {
struct SpinDensityMatrix {
    Eigen::Vector3d bp_matrix, bm_matrix;
    Eigen::Matrix3d c_matrix;
    double norm_factor;
};

SpinDensityMatrix spinDensityMatrix(double sqrt_s_hat, double cos_th);
}

#endif  // SRC_SPIN_DENSITY_H
