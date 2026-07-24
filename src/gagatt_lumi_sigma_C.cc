// gagatt_lumi_sigma_C.cc
//
// Standalone tool to compute <C[rho]> and sigma_C as a function of
// integrated luminosity L, for the significance plot of concurrence.
//
// Output columns: L [ab^-1],  <C[rho]>,  sigma_C

#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <string>
#include <vector>
#include "constants.h"
#include "mc.h"
#include "mc_helper.h"
#include "photon.h"
#include "spin_density.h"

using namespace gagatt;

void printUsage(const char *prog) {
    std::cerr << "usage: " << prog
              << " <output> <sqrt_s> <pe1> <pe2>\n"
                 "   [n_events] [L_scan_min_ab] [L_scan_max_ab] [n_L_points]\n"
                 "   [--coscut <val>] [--sqrts_max <val>]\n\n"
              << " output       : output file name\n"
              << " sqrt_s       : e+e- CM energy [GeV] (e.g. 500 or 1000)\n"
              << " pe1          : electron polarization\n"
              << " pe2          : positron polarization\n"
              << " n_events     : MC events per seed [default 10_000_000]\n"
              << " L_scan_min   : scan start [ab^-1] [default 0.001]\n"
              << " L_scan_max   : scan end [ab^-1] [default 0.5]\n"
              << " n_L_points   : number of luminosity pts [default 100]\n"
              << " --coscut     : restrict |cos(Theta)| < val\n"
              << " --sqrts_max  : upper edge of sqrt(s_hat) scan [GeV]\n";
}

int main(int argc, char *argv[]) {
    // Parse --coscut and --sqrts_max flags first (same pattern as gagatt_mc.cc)
    MCConfig cfg;

    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);

        auto nextVal = [&](const char *flag) -> const char * {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: " << flag << " requires a value\n";
                std::exit(EXIT_FAILURE);
            }
            return argv[++i];
        };

        if (arg == "--coscut") {
            const double v = std::stod(nextVal("--coscut"));
            if (v <= 0.0 || v >= 1.0) {
                std::cerr << std::format(
                    "ERROR: --coscut must be in (0, 1); got {:.4f}\n", v);
                return EXIT_FAILURE;
            }
            cfg.cos_th_min = -v;
            cfg.cos_th_max = v;
        } else if (arg == "--sqrts_max") {
            cfg.sqrts_max = std::stod(nextVal("--sqrts_max"));
        } else {
            pos.push_back(argv[i]);
        }
    }

    // positional args: output, sqrt_s, pe1, pe2, [n_events, ...]
    if (pos.size() < 4 || pos.size() > 8) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const std::string output_file = pos[0];
    cfg.sqrt_s = std::stod(pos[1]);
    cfg.pe1 = std::stod(pos[2]);
    cfg.pe2 = std::stod(pos[3]);
    if (pos.size() >= 5) { cfg.n_events = std::stoll(pos[4]); }
    if (pos.size() >= 6) { cfg.L_scan_min_ab = std::stod(pos[5]); }
    if (pos.size() >= 7) { cfg.L_scan_max_ab = std::stod(pos[6]); }
    if (pos.size() >= 8) { cfg.n_L_points = std::stoi(pos[7]); }

    // sqrts_max defaults to sqrt_s unless --sqrts_max was given.
    if (cfg.sqrts_max < 0.0) { cfg.sqrts_max = cfg.sqrt_s; }
    if (cfg.sqrts_max > cfg.sqrt_s) {
        std::cerr << std::format(
            "ERROR: --sqrts_max ({:.1f}) cannot exceed sqrt_s ({:.1f})\n",
            cfg.sqrts_max, cfg.sqrt_s);
        return EXIT_FAILURE;
    }

    // --- Build weight table (same as runMC) ---
    const double pc1 = -cfg.pe1;
    const double pc2 = -cfg.pe2;
    const double sqrts_min = cfg.sqrts_min;
    const double d_sqrts = (cfg.sqrts_max - sqrts_min) / cfg.n_sqrts;
    const double d_cos = (cfg.cos_th_max - cfg.cos_th_min) / cfg.n_cos;

    auto zcache = buildLumiCache(cfg, pc1, pc2, sqrts_min, d_sqrts);
    auto wt = buildWeightTable(cfg, zcache, sqrts_min, d_sqrts, d_cos);

    if (wt.bin_weights.empty() || wt.total_weight <= 0.0) {
        std::cerr << "ERROR: weight table is empty — check kinematics.\n";
        return EXIT_FAILURE;
    }

    // sigma_eff = sigma * BR(t tbar -> l+ l-)
    const double sigma_eff_fb = wt.total_weight * BRLL;

    std::cout << std::format(
        "-- sigma_tot = {:.6f} fb, sigma_eff = {:.6f} fb\n", wt.total_weight,
        sigma_eff_fb);
    std::cout << std::format(
        "-- theory C = {:.6f}, theory D = {:.6f}, theory m12 = {:.6f}\n",
        wt.theory_concurrence, wt.theory_D, wt.theory_m12);

    // --- Luminosity scan ---
    const int K = std::max(1, cfg.n_lumi_seeds);
    const double dL = (cfg.L_scan_max_ab - cfg.L_scan_min_ab) /
                      static_cast<double>(cfg.n_L_points - 1);

    std::mt19937_64 rng(cfg.seed ? cfg.seed : std::random_device{}());

    std::cout << std::format(
        "\n-- luminosity scan [{:.3f}, {:.3f}] ab^-1, {} points, {} "
        "seeds/point\n",
        cfg.L_scan_min_ab, cfg.L_scan_max_ab, cfg.n_L_points, K);
    std::cout << std::format(" {:>6} {:>10} {:>12} {:>12} {:>12}\n", "pt",
                             "L[ab^-1]", "N(L)", "mc_C", "sigma_C");

    std::ofstream ofs(output_file);
    if (!ofs) {
        std::cerr << "ERROR: cannot open " << output_file << " for writing\n";
        return EXIT_FAILURE;
    }

    // Header
    ofs << std::format("# L, mc_C, sigma_C vs luminosity\n");
    ofs << std::format("# sqrt_s = {:.1f} GeV\n", cfg.sqrt_s);
    ofs << std::format("# pe1 = {:+.2f}\n", cfg.pe1);
    ofs << std::format("# pe2 = {:+.2f}\n", cfg.pe2);
    if (cfg.cos_th_min != -1.0 || cfg.cos_th_max != 1.0) {
        ofs << std::format("# cos_th_min = {:.4f}\n", cfg.cos_th_min);
        ofs << std::format("# cos_th_max = {:.4f}\n", cfg.cos_th_max);
    }
    if (cfg.sqrts_max != cfg.sqrt_s) {
        ofs << std::format("# sqrts_max = {:.1f} GeV\n", cfg.sqrts_max);
    }
    ofs << std::format("# n_lumi_seeds = {}\n", K);
    ofs << std::format("# sigma_tot = {:.6f} fb\n", wt.total_weight);
    ofs << std::format("# sigma_eff = {:.6f} fb\n", sigma_eff_fb);
    ofs << std::format("# theory C = {:.6f}\n", wt.theory_concurrence);
    ofs << std::format("#\n");
    ofs << std::format("# {:>12s} {:>16s} {:>16s}\n", "L[ab^-1]", "mc_C",
                       "sigma_C");

    for (int p = 0; p < cfg.n_L_points; ++p) {
        const double L_ab = cfg.L_scan_min_ab + p * dL;
        const double L_fb = L_ab * 1.0e3;        // ab^-1 -> fb^-1
        const double N_L = sigma_eff_fb * L_fb;  // physical event count
        const long long N_events = static_cast<long long>(std::max(1.0, N_L));

        if (N_L <= 0.0) {
            ofs << std::format(" {:12.6f} {:16s} {:16s}\n", L_ab, "N/A", "N/A");
            std::cout << std::format(
                " {:>6} {:>10.4f} {:>12.1f} {:>12} {:>12} {:>12}\n", p + 1,
                L_ab, 0.0, "N/A", "N/A", "N/A");
            continue;
        }

        std::vector<double> mc_C_vec(K), sigma_C_vec(K);
        for (int k = 0; k < K; ++k) {
            const uint64_t seed_k = rng();
            std::mt19937_64 rng_k(seed_k);

            const EventLoopResult ev_k =
                runEventLoop(N_events, wt, rng_k, /*verbose=*/false);
            const ReconstructedMC recon_k = reconstructFromMoments(ev_k);

            mc_C_vec[k] = recon_k.mc_concurrence;
            sigma_C_vec[k] = recon_k.sigma_concurrence;
        }

        const double mean_mc_C =
            std::accumulate(mc_C_vec.begin(), mc_C_vec.end(), 0.0) / K;
        const double mean_sigma_C =
            std::accumulate(sigma_C_vec.begin(), sigma_C_vec.end(), 0.0) / K;

        ofs << std::format(" {:12.6f} {:16.8e} {:16.8e}\n", L_ab, mean_mc_C,
                           mean_sigma_C);

        std::cout << std::format(
            " {:>6} {:>10.4f} {:>12.1f} {:>12.6f} {:>12.6f}\n", p + 1, L_ab,
            N_L, mean_mc_C, mean_sigma_C);
    }

    ofs.close();
    std::cout << std::format("\n-- output written to: {}\n", output_file);

    return EXIT_SUCCESS;
}
