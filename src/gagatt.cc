#include <cmath>
#include <iostream>
#include <numbers>

#include "amplitude.h"
#include "constants.h"

using H = gagatt::Helicity;

int main() {
    std::cout << "Hello, world!\n";
    std::cout << "M(top) = " << gagatt::MTOP << '\n';

    const double sqrt_s_hat = 400.0;
    const double cos_th = std::cos(45.0 * std::numbers::pi / 180.0);
    std::cout << "sqrt(s) = " << sqrt_s_hat << '\n';
    std::cout << "cos(theta) = " << cos_th << '\n';

    auto coeffs = gagatt::computeCoeffsOnShell(sqrt_s_hat, cos_th);

    std::cout << "C1(on-shell) = " << coeffs.c1 << '\n';
    std::cout << "C2(on-shell) = " << coeffs.c2 << '\n';
    std::cout << "C3(on-shell) = " << coeffs.c3 << '\n';
    std::cout << "C4(on-shell) = " << coeffs.c4 << '\n';
}
