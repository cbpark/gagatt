#include <cmath>
#include <iostream>
#include <numbers>

#include "amplitude.h"
#include "constants.h"
#include "spin_density.h"

using H = gagatt::Helicity;

int main() {
    std::cout << "Hello, world!\n";
    std::cout << "M(top) = " << gagatt::MTOP << '\n';

    const double sqrt_s_hat = 400.0;
    const double cos_th = std::cos(45.0 * std::numbers::pi / 180.0);
    std::cout << "sqrt(s) = " << sqrt_s_hat << '\n';
    std::cout << "cos(theta) = " << cos_th << '\n';

    auto coeffs = gagatt::computePolCoeffs(sqrt_s_hat, cos_th);

    std::cout << "C1(on-shell) = " << coeffs.c1 << '\n';
    std::cout << "C2(on-shell) = " << coeffs.c2 << '\n';
    std::cout << "C3(on-shell) = " << coeffs.c3 << '\n';
    std::cout << "C4(on-shell) = " << coeffs.c4 << '\n';
    std::cout << "C5(on-shell) = " << coeffs.c5 << '\n';
    std::cout << "C6(on-shell) = " << coeffs.c6 << '\n';
    std::cout << "C7(on-shell) = " << coeffs.c7 << '\n';
    std::cout << "C8(on-shell) = " << coeffs.c8 << '\n';
    std::cout << "C9(on-shell) = " << coeffs.c9 << '\n';
    std::cout << "C10(on-shell) = " << coeffs.c10 << '\n';
    std::cout << "C11(on-shell) = " << coeffs.c11 << '\n';
    std::cout << "C12(on-shell) = " << coeffs.c12 << '\n';
    std::cout << "C13(on-shell) = " << coeffs.c13 << '\n';
    std::cout << "C14(on-shell) = " << coeffs.c14 << '\n';
    std::cout << "C15(on-shell) = " << coeffs.c15 << '\n';
    std::cout << "C16(on-shell) = " << coeffs.c16 << '\n';

    // std::cout << gagatt::I2I2 << '\n';

    auto sdc = gagatt::SDMatrixCoefficients(sqrt_s_hat, cos_th);
    std::cout << "A: " << sdc.norm_factor << '\n';
    std::cout << "C(ij):\n" << sdc.cc << '\n';

    auto rho = gagatt::spinDensityMatrix(sqrt_s_hat, cos_th);
    std::cout << "rho:\n " << rho << '\n';
}
