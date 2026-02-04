#include <array>
#include <fstream>
#include <iostream>
#include "photon.h"

int main() {
    const double x = 4.8;

    std::array<double, 100> z{}, lumi{};
    std::ofstream fout("tmp/lumi.dat");
    for (int i = 0; i < 100; ++i) {
        z[i] = i / 100.0;
        lumi[i] = gagatt::photonLuminosityUnpol(z[i], x, 1.0, -1.0, 1.0, -1.0);
        // lumi[i] = gagatt::photonLuminosityUnpol(z[i], x, 0, 0, 0, 0);

        fout << z[i] << "  " << lumi[i] << '\n';
    }
}
