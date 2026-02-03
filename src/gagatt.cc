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

    auto on_shell_amp = gagatt::onShellAmp(sqrt_s_hat, cos_th, H::PLUS, H::PLUS,
                                           H::PLUS, H::MINUS);
    std::cout << "on_shell_amp = " << on_shell_amp << '\n';

    const double c1_on_shell = gagatt::c1OnShell(sqrt_s_hat, cos_th);
    std::cout << "C1(on-shell) = " << c1_on_shell << '\n';

    const double c3_on_shell = gagatt::c3OnShell(sqrt_s_hat, cos_th);
    std::cout << "C3(on-shell) = " << c3_on_shell << '\n';
}
