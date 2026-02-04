#include <array>
#include <fstream>
#include <iostream>
#include "photon.h"

int main() {
    const double x = 4.8;

    std::array<double, 100> y{}, f{};
    std::ofstream fout("f.dat");
    for (int i = 0; i < 100; ++i) {
        y[i] = i / 100.0;
        f[i] = gagatt::fLumi(x, y[i], 1.0, -1.0);

        fout << y[i] << "  " << f[i] << '\n';
    }
}
