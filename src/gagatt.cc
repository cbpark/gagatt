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

    std::cout << "C1: " << coeffs.c1 << "$\n";
    std::cout << "C2: " << coeffs.c2 << "$\n";
    std::cout << "C3: " << coeffs.c3 << "$\n";
    std::cout << "C4: " << coeffs.c4 << "$\n";
    std::cout << "C5: " << coeffs.c5 << "$\n";
    std::cout << "C6: " << coeffs.c6 << "$\n";
    std::cout << "C7: " << coeffs.c7 << "$\n";
    std::cout << "C8: " << coeffs.c8 << "$\n";
    std::cout << "C9: " << coeffs.c9 << "$\n";
    std::cout << "C10: " << coeffs.c10 << "$\n";
    std::cout << "C11: " << coeffs.c11 << "$\n";
    std::cout << "C12: " << coeffs.c12 << "$\n";
    std::cout << "C13: " << coeffs.c13 << "$\n";
    std::cout << "C14: " << coeffs.c14 << "$\n";
    std::cout << "C15: " << coeffs.c15 << "$\n";
    std::cout << "C16: " << coeffs.c16 << "$\n";

    auto sdc = gagatt::SDMatrixCoefficients(sqrt_s_hat, cos_th);
    // std::cout << "A: " << sdc.norm_factor << '\n';
    // std::cout << "C(nn): " << sdc.c_nn() << '\n';
    // std::cout << "C(ij):\n" << sdc.cc << '\n';

    auto rho = gagatt::spinDensityMatrix(sdc);
    std::cout << "rho:\n " << rho << '\n';

    auto rho_t2 = gagatt::partialTransposeB(rho);
    std::cout << "rho(T2):\n " << rho_t2 << '\n';

    std::cout << "PH criterion: " << gagatt::isEntangledByPH(rho) << '\n';

    std::cout << "Concurrence: " << gagatt::getConcurrence(rho) << '\n';
    std::cout << "Conc criterion: " << gagatt::isEntangledByConcurrence(rho)
              << '\n';

    std::cout << "Violation of Bell inequality: "
              << gagatt::violatesBellInequality(sdc) << '\n';

    std::cout << "Tr[C]: " << gagatt::isEntangledByD(sdc) << '\n';
}
