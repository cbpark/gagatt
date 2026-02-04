#include <array>
#include <fstream>
#include <iostream>
#include "amplitude.h"
#include "photon.h"

int main() {
    const double x = 4.8;

    std::array<double, 100> z{}, lumi{}, lumi_p{}, lumi_m;
    std::ofstream fout("data/photon_luminosity.dat");
    for (int i = 0; i < 100; ++i) {
        z[i] = i / 100.0;

        lumi[i] = gagatt::photonLuminosityUnpol(z[i], x, 1.0, -1.0, 1.0, -1.0);
        lumi_p[i] = gagatt::photonLuminosity(z[i], x, 1.0, -1.0, 1.0, -1.0,
                                             {gagatt::Helicity::PLUS},
                                             {gagatt::Helicity::PLUS});
        lumi_m[i] = gagatt::photonLuminosity(z[i], x, -1.0, 1.0, -1.0, 1.0,
                                             {gagatt::Helicity::PLUS},
                                             {gagatt::Helicity::PLUS});

        fout << z[i] << "  " << lumi[i] << "  " << lumi_p[i] << "  "
             << lumi_m[i] << '\n';
    }
}
