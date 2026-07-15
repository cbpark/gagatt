#ifndef SRC_MC_HELPER_H
#define SRC_MC_HELPER_H

#include <random>
#include <utility>
#include <vector>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include "mc.h"
#include "photon.h"
#include "spin_density.h"

namespace gagatt {
// -----------------------------------------------------------------------
// ZCacheEntry / LumiCache
//
// Precomputed (LumiWeights, L_tot) per sqrt_s_hat bin, indexed by j in
// [0, n_sqrts). Reused by buildWeightTable() so lumiWeightsAndTotal() is
// evaluated exactly once per sqrt_s_hat bin instead of once per (i, j) cell.
// -----------------------------------------------------------------------
struct ZCacheEntry {
    LumiWeights lw;
    double L_tot;
};

// Builds the z-cache over n_sqrts bins spanning [sqrts_min, sqrts_max].
// pc1 = -pe1, pc2 = -pe2 (PePc = -1 convention, see runMC).
std::vector<ZCacheEntry> buildLumiCache(const MCConfig &cfg, double pc1,
                                        double pc2, double sqrts_min,
                                        double d_sqrts);

// -----------------------------------------------------------------------
// partialXsec / eventRate
//
// d sigma_hat / d cos_th (helicity-summed, luminosity-weighted) and the
// differential event rate d^2 sigma / (d sqrt_s_hat d cos_th) [fb/GeV].
// See mc_helper.cc for the full derivation/unit bookkeeping.
// -----------------------------------------------------------------------
double partialXsec(double sqrt_s_hat, double cos_th,
                   const SDMatrixCoefficients &sdc);

double eventRate(double sqrt_s_hat, double cos_th,
                 const SDMatrixCoefficients &sdc, double L_tot, double sqrt_s);

// -----------------------------------------------------------------------
// WeightTable
//
// Flattened 2-D (cos_th, sqrt_s_hat) weight table plus the cached
// SDMatrixCoefficients for each bin (so the event loop doesn't need to
// recompute sdc for the sampled bin), and the luminosity+phase-space
// weighted theory averages of D, Tr[C], negativity, concurrence, m12.
// -----------------------------------------------------------------------
struct WeightTable {
    std::vector<double> bin_weights;
    std::vector<SDMatrixCoefficients> sdc_cache;
    double total_weight = 0.0;

    double theory_tr_c = 0.0;
    double theory_D = 0.0;
    double theory_negativity = 0.0;
    double theory_concurrence = 0.0;
    double theory_m12 = 0.0;
};

// Builds the flattened weight table over the (n_cos x n_sqrts) grid and
// accumulates the theory averages in a single pass. Returns an empty
// bin_weights vector (size 0) if total_weight <= 0 (caller must check).
WeightTable buildWeightTable(const MCConfig &cfg,
                             const std::vector<ZCacheEntry> &zcache,
                             double sqrts_min, double d_sqrts, double d_cos);

// -----------------------------------------------------------------------
// sampleDecayAngles
//
// Draw (q+, q-) unit vectors for l+/l- in the top/anti-top rest frames
// from the joint angular distribution via accept/reject sampling.
// See mc_helper.cc for the distribution and envelope details.
// -----------------------------------------------------------------------
std::pair<Eigen::Vector3d, Eigen::Vector3d> sampleDecayAngles(
    const SDMatrixCoefficients &sdc, std::mt19937_64 &rng);

// -----------------------------------------------------------------------
// EventLoopResult
//
// Raw first/second moments of the outer product q+_i * q-_j, accumulated
// over n_accepted sampled events. Kept unprocessed (not yet turned into
// C_ij) so the same accumulators can be reused for binned (D_vs_sqrts,
// D_vs_costh) scans in the future.
// -----------------------------------------------------------------------
struct EventLoopResult {
    Eigen::Matrix3d S1_qpqm = Eigen::Matrix3d::Zero();  // sum q+_i q-_j
    Eigen::Matrix3d S2_qpqm = Eigen::Matrix3d::Zero();  // sum (q+_i q-_j)^2
    long long n_accepted = 0;
};

// Runs the MC event loop: draws a (cos_th, sqrt_s_hat) bin proportional to
// wt.bin_weights, samples decay angles, and accumulates moments.
// Prints progress to stdout every ~10% of cfg.n_events.
EventLoopResult runEventLoop(const MCConfig &cfg, const WeightTable &wt,
                             std::mt19937_64 &rng);

// -----------------------------------------------------------------------
// ReconstructedMC
//
// All MC-side quantities derived from the raw event-loop moments:
// C_ij, its uncertainty, D, Tr[C], reconstructed rho and its entanglement
// measures, and the Horodecki m12 parameter with its (rough) uncertainty.
// -----------------------------------------------------------------------
struct ReconstructedMC {
    Eigen::Matrix3d mc_cij = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d sigma_cij = Eigen::Matrix3d::Zero();

    double mc_tr_c = 0.0;
    double sigma_tr_c = 0.0;

    double mc_D = 0.0;
    double sigma_D = 0.0;
    double significance_D = 0.0;  // (-1/3 - D) / sigma_D when D < -1/3

    double mc_negativity = 0.0;
    double mc_concurrence = 0.0;

    double mc_m12 = 0.0;
    double sigma_m12 = 0.0;
    double significance_bell = 0.0;  // (m12 - 1) / sigma_m12 when m12 > 1
};

// Turns raw event-loop moments into C_ij and all derived quantities.
ReconstructedMC reconstructFromMoments(const EventLoopResult &ev);

// -----------------------------------------------------------------------
// Luminosity scan
//
// Rescales significance_D / significance_bell (computed at N_MC events)
// to an arbitrary integrated luminosity L via statistical sqrt(N) scaling.
// See mc.cc / mc_helper.cc for the full derivation.
// -----------------------------------------------------------------------
std::vector<LumiScanPoint> computeLumiScan(const MCConfig &cfg,
                                           double sigma_eff_fb,
                                           long long n_accepted,
                                           double significance_D, double mc_m12,
                                           double sigma_m12);
}  // namespace gagatt

#endif  // SRC_MC_HELPER_H
