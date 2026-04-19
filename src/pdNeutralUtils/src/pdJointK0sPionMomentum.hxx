#ifndef pdJointK0sPionMomentum_hxx
#define pdJointK0sPionMomentum_hxx

#include "pdDataClasses.hxx"
#include "TVector3.h"
#include <limits>
#include <vector>

struct JointK0sPionMomentumGridResult {
  bool ok = false;
  Float_t p1 = -999.f;
  Float_t p2 = -999.f;
  Float_t bestScore = -999.f;
  Float_t invMassAtBest = -999.f;
  /// Post-fit only: |C_mass| / max(|logL1+logL2|, eps); C_mass = −(M−m_K0)²/(2σ²) (mass term without penaltyScale).
  /// Large R ⇒ mass constraint dominated the optimum vs track TLE likelihoods.
  Float_t constraintRatioR = -999.f;
};

/// Filled by ComputeSigmaMEventGeV for debugging / tree export (GeV; NaN if not computed).
struct JointK0sSigmaMEventDiagnostics {
  double sigma_p1_gev = std::numeric_limits<double>::quiet_NaN();
  double sigma_p2_gev = std::numeric_limits<double>::quiet_NaN();
  double sigma_m_event_gev = std::numeric_limits<double>::quiet_NaN();
  double dm_dp1 = std::numeric_limits<double>::quiet_NaN();
  double dm_dp2 = std::numeric_limits<double>::quiet_NaN();
  double cos_theta = std::numeric_limits<double>::quiet_NaN();
};

namespace pdJointK0sPionMomentum {

/// Piecewise-linear log L(p) with clamp to endpoints (same scan as free-range ML).
double InterpolateLogLikelihoodClamped(const std::vector<double>& pAxis, const std::vector<double>& logL, double p);

/// Invariant mass of two pions with magnitudes a1,a2 and directions (any length; unit used internally).
bool PionPairInvariantMassGeV(double a1, double a2, const TVector3& dir1, const TVector3& dir2, double& massGeV);

/// π+π− invariant mass and ∂m/∂p1, ∂m/∂p2 for fixed cosθ = û1·û2 (p1,p2 magnitudes in GeV/c).
bool PionPiPiInvariantMassAndDerivativesGeV(double p1, double p2, double cos_theta, double mpi_gev, double& m_gev,
                                           double& dm_dp1, double& dm_dp2);

/// σ_p from discrete TLE log-likelihood curve: NLL = -log L; curvature finite-difference, else ΔNLL=0.5 width.
double EstimateSigmaPFromLogLikelihoodCurve(const std::vector<double>& pAxis, const std::vector<double>& logL,
                                             double p0_gev);

/// Event σ_m from error propagation at marginal TLE maxima (p1_0,p2_0); floor/cap/inflate near-collinear.
/// fallback_sigma_m_gev used if propagation is non-finite (returns true anyway with safe value).
bool ComputeSigmaMEventGeV(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                           const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                           const TVector3& dir1, const TVector3& dir2, double fallback_sigma_m_gev,
                           JointK0sSigmaMEventDiagnostics& out);

/// Maximize logL1(p1)+logL2(p2) - 0.5*penaltyScale*((M-mK0s)/sigmaMassGeV)^2 on a coarse grid, optional refinement.
JointK0sPionMomentumGridResult FitJointMomentaOnGrid(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                                                     const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                                                     const TVector3& dir1, const TVector3& dir2, double pMinGeV,
                                                     double pMaxGeV, double pStepGeV, double mK0sGeV, double sigmaMassGeV,
                                                     double penaltyScale, int refineFactor);

} // namespace pdJointK0sPionMomentum

#endif
