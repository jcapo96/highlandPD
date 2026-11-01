#ifndef pdDedxPdfModel_h
#define pdDedxPdfModel_h

// Likelihood model of the measured dE/dx of a hit, and the TLE (free-range) scan built on it.
//
// Self-contained (ROOT only), so that the same code runs in highlandPD and from PyROOT on the pion dumps
// (claude_work/scripts/k0s_pionmom/py). Pion-momentum audit, plans/pion_momentum_audit.md, items B1-B3, B8, B10:
//   - the per-hit pitch enters the Landau-Vavilov width and most probable value (B3);
//   - the energy-loss distribution is convolved with a Gaussian detector resolution (B2);
//   - the PDF is renormalised to the dE/dx window used to select the hits, and a uniform outlier component keeps
//     every hit in the likelihood, so no hit is skipped for having zero probability (B1);
//   - the Gaussian (kappa > 10) branch is normalised (B8);
//   - Bethe-Bloch, the CSDA range table and the PDF use one liquid-argon density (B10), the same as the Python
//     range momenta (lambda_feasibility.p_from_range).
// The model is off by default (Model::enabled = false); the callers then keep their legacy code path unchanged.

#include <vector>

namespace pdDedxPdf {

struct Model {
  bool enabled = false;
  bool perHitPitch = true;     // use the hit pitch; fixedPitchCm when false or when the hit pitch is invalid
  bool trackMedianPitch = false;  // with perHitPitch: the median pitch of the selected hits for every hit (ScanTLE)
  double fixedPitchCm = 0.65;
  double resolutionFrac = 0.0; // Gaussian sigma as a fraction of the Bethe-Bloch mean dE/dx
  double resolutionAbs = 0.0;  // MeV/cm, added in quadrature
  double outlierFrac = 0.0;    // weight of a uniform component over the dE/dx window
  bool truncate = true;        // renormalise the PDF to [dedxMin, dedxMax]
  double dedxMin = 0.5;        // MeV/cm, the hit-selection window of the TLE
  double dedxMax = 5.0;
  double densityGcm3 = 1.3954;
  double dedxScale = 1.0;      // multiplies the measured dE/dx in ScanTLE (calibration of the data to the model)
};

void SetModel(const Model& m);  // also clears the PDF cache
const Model& GetModel();
bool Enabled();

// Bethe-Bloch mean dE/dx in liquid argon with the Sternheimer density effect [MeV/cm].
double MeanDedx(double keMeV, double massMeV);
// CSDA range from Bethe-Bloch (0.1 MeV - 5 GeV table, linear extrapolation above) and its inverse.
double CsdaRangeCm(double keMeV, double massMeV);
double KeFromRangeCm(double rangeCm, double massMeV);
double MomentumGeVFromRangeCm(double rangeCm, double massMeV);

// ln of the PDF of the measured dE/dx of one hit, for a particle of kinetic energy keMeV crossing pitchCm.
// Always finite with the model enabled (floored at ln 1e-12 if the outlier fraction is zero).
double LogPdf(double dedx, double keMeV, double massMeV, double pitchCm);

// TLE scan: hypothesis "the particle stops L cm after the hit with residual range 0", L from -min(rr) to Lmax.
// Hits are the already selected ones (skips and dE/dx window applied). Momentum of each hypothesis from the CSDA
// range trackLengthCm + L. Refines with fineStepCm when any hypothesis is below lowPRefineGeV (as the legacy scan).
struct ScanResult {
  std::vector<double> offsetCm, pGeV, logL;
  int best = -1;
};
bool ScanTLE(const std::vector<double>& rr, const std::vector<double>& dedx, const std::vector<double>& pitch,
             double trackLengthCm, double massMeV, double LmaxCm, double stepCm, double fineStepCm,
             double lowPRefineGeV, ScanResult& out);

}  // namespace pdDedxPdf

#endif
