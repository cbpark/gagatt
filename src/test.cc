#include <Eigen/Dense>
#include <iostream>

int main() {
    Eigen::Matrix2d m = Eigen::Matrix2d::Identity();
    std::cout << "Eigen is working! Matrix:\n" << m << std::endl;
}
