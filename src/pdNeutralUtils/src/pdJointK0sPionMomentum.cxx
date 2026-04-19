#include "pdJointK0sPionMomentum.hxx"
#include <algorithm>
#include <cmath>
#include <limits>

namespace pdJointK0sPionMomentum {
namespace {

constexpr double kFloorSigmaMGeV = 0.005;  // 5 MeV
constexpr double kCapSigmaMGeV = 0.05;    // 50 MeV
constexpr double kCollinearCosCut = 0.95;
constexpr double kCollinearInflate = 2.0;

double NLLFromLogLikelihoodCurve(const std::vector<double>& pAxis, const std::vector<double>& logL, double p) {
  const double ll = InterpolateLogLikelihoodClamped(pAxis, logL, p);
  if (!std::isfinite(ll)) return std::numeric_limits<double>::infinity();
  return -ll;
}

double MedianPositiveSpacing(const std::vector<double>& p) {
  if (p.size() < 2u) return 1e-4;
  std::vector<double> d;
  d.reserve(p.size() - 1u);
  for (size_t i = 1; i < p.size(); ++i) {
    const double s = p[i] - p[i - 1];
    if (s > 0. && std::isfinite(s)) d.push_back(s);
  }
  if (d.empty()) return 1e-4;
  const size_t mid = d.size() / 2u;
  std::nth_element(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(mid), d.end());
  return std::max(1e-9, d[mid]);
}

/// One-sided distance from p0 (in momentum) where NLL first reaches N0 + 0.5 (linear scan then refine).
double OneSidedHalfDeltaNLLWidth(const std::vector<double>& pAxis, const std::vector<double>& logL, double p0,
                                 int direction /* +1 right, -1 left */) {
  const double N0 = NLLFromLogLikelihoodCurve(pAxis, logL, p0);
  if (!std::isfinite(N0)) return -1.;
  const double target = N0 + 0.5;
  const double pLo = pAxis.front();
  const double pHi = pAxis.back();
  const int nScan = 400;
  if (direction > 0) {
    for (int k = 1; k <= nScan; ++k) {
      const double p = p0 + (pHi - p0) * (static_cast<double>(k) / static_cast<double>(nScan));
      if (p > pHi) break;
      const double Np = NLLFromLogLikelihoodCurve(pAxis, logL, p);
      if (std::isfinite(Np) && Np >= target) return p - p0;
    }
  } else {
    for (int k = 1; k <= nScan; ++k) {
      const double p = p0 - (p0 - pLo) * (static_cast<double>(k) / static_cast<double>(nScan));
      if (p < pLo) break;
      const double Np = NLLFromLogLikelihoodCurve(pAxis, logL, p);
      if (std::isfinite(Np) && Np >= target) return p0 - p;
    }
  }
  return -1.;
}

double SigmaFromDeltaNLLHalf(const std::vector<double>& pAxis, const std::vector<double>& logL, double p0,
                             double fallbackWidth) {
  const double wR = OneSidedHalfDeltaNLLWidth(pAxis, logL, p0, +1);
  const double wL = OneSidedHalfDeltaNLLWidth(pAxis, logL, p0, -1);
  if (wR > 0. && wL > 0.) return 0.5 * (wR + wL);
  if (wR > 0.) return wR;
  if (wL > 0.) return wL;
  return fallbackWidth;
}

double MassPenaltyChi2(double massGeV, double mK0sGeV, double sigmaMassGeV, double penaltyScale) {
  if (!(penaltyScale > 0.) || !std::isfinite(penaltyScale)) return 0.;
  const double sig = (sigmaMassGeV > 1e-9 && std::isfinite(sigmaMassGeV)) ? sigmaMassGeV : 1e-3;
  if (!std::isfinite(massGeV) || !std::isfinite(mK0sGeV)) return 1.e9;
  const double d = (massGeV - mK0sGeV) / sig;
  return 0.5 * penaltyScale * d * d;
}

void ScanRectangle(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                   const std::vector<double>& p2Axis, const std::vector<double>& logL2, const TVector3& dir1,
                   const TVector3& dir2, double p1Lo, double p1Hi, double p2Lo, double p2Hi, double pStepGeV,
                   double mK0sGeV, double sigmaMassGeV, double penaltyScale, double& bestS, double& bestP1,
                   double& bestP2, double& bestM) {
  if (!(pStepGeV > 0.) || !(p1Hi >= p1Lo) || !(p2Hi >= p2Lo)) return;
  for (double a1 = p1Lo; a1 <= p1Hi + 1e-12; a1 += pStepGeV) {
    const double L1 = InterpolateLogLikelihoodClamped(p1Axis, logL1, a1);
    if (!std::isfinite(L1)) continue;
    for (double a2 = p2Lo; a2 <= p2Hi + 1e-12; a2 += pStepGeV) {
      const double L2 = InterpolateLogLikelihoodClamped(p2Axis, logL2, a2);
      if (!std::isfinite(L2)) continue;
      double mPiPi = 0.;
      if (!PionPairInvariantMassGeV(a1, a2, dir1, dir2, mPiPi)) continue;
      const double pen = MassPenaltyChi2(mPiPi, mK0sGeV, sigmaMassGeV, penaltyScale);
      const double S = L1 + L2 - pen;
      if (std::isfinite(S) && S > bestS) {
        bestS = S;
        bestP1 = a1;
        bestP2 = a2;
        bestM = mPiPi;
      }
    }
  }
}

} // namespace

double InterpolateLogLikelihoodClamped(const std::vector<double>& pAxis, const std::vector<double>& logL, double p) {
  if (pAxis.size() < 2u || pAxis.size() != logL.size()) return std::numeric_limits<double>::quiet_NaN();
  if (!std::isfinite(p)) return std::numeric_limits<double>::quiet_NaN();
  const double pLo = pAxis.front();
  const double pHi = pAxis.back();
  if (p <= pLo) return logL.front();
  if (p >= pHi) return logL.back();

  auto it = std::upper_bound(pAxis.begin(), pAxis.end(), p);
  const size_t j = static_cast<size_t>(std::distance(pAxis.begin(), it));
  const size_t i = j - 1u;
  const double p0 = pAxis[i];
  const double p1 = pAxis[j];
  const double t = (p - p0) / (p1 - p0);
  return (1. - t) * logL[i] + t * logL[j];
}

bool PionPairInvariantMassGeV(double a1, double a2, const TVector3& dir1, const TVector3& dir2, double& massGeV) {
  constexpr double kPionMassGeV = 0.13957;
  if (!(a1 > 0.) || !(a2 > 0.) || !std::isfinite(a1) || !std::isfinite(a2)) return false;
  if (dir1.Mag2() <= 0. || dir2.Mag2() <= 0.) return false;
  const TVector3 u1 = dir1.Unit();
  const TVector3 u2 = dir2.Unit();
  const TVector3 p3v = a1 * u1 + a2 * u2;
  const double e1 = std::sqrt(a1 * a1 + kPionMassGeV * kPionMassGeV);
  const double e2 = std::sqrt(a2 * a2 + kPionMassGeV * kPionMassGeV);
  const double eTot = e1 + e2;
  const double m2 = eTot * eTot - p3v.Mag2();
  if (!(m2 > 0.) || !std::isfinite(m2)) return false;
  massGeV = std::sqrt(m2);
  return std::isfinite(massGeV) && massGeV > 0.;
}

bool PionPiPiInvariantMassAndDerivativesGeV(double p1, double p2, double cos_theta, double mpi_gev, double& m_gev,
                                           double& dm_dp1, double& dm_dp2) {
  if (!(p1 > 0.) || !(p2 > 0.) || !std::isfinite(p1) || !std::isfinite(p2) || !std::isfinite(mpi_gev)) return false;
  const double ct = std::max(-1.0, std::min(1.0, cos_theta));
  const double e1 = std::sqrt(p1 * p1 + mpi_gev * mpi_gev);
  const double e2 = std::sqrt(p2 * p2 + mpi_gev * mpi_gev);
  if (!(e1 > 0.) || !(e2 > 0.)) return false;
  const double m2 = 2.0 * mpi_gev * mpi_gev + 2.0 * (e1 * e2 - p1 * p2 * ct);
  if (!(m2 > 0.) || !std::isfinite(m2)) return false;
  m_gev = std::sqrt(m2);
  if (!(m_gev > 0.) || !std::isfinite(m_gev)) return false;
  dm_dp1 = ((p1 * e2 / e1) - (p2 * ct)) / m_gev;
  dm_dp2 = ((p2 * e1 / e2) - (p1 * ct)) / m_gev;
  return std::isfinite(dm_dp1) && std::isfinite(dm_dp2);
}

double EstimateSigmaPFromLogLikelihoodCurve(const std::vector<double>& pAxis, const std::vector<double>& logL,
                                             double p0_gev) {
  const double medianStep = MedianPositiveSpacing(pAxis);
  const double span = pAxis.back() - pAxis.front();
  const double fallbackW = std::max(medianStep * 3.0, std::max(0.05 * std::max(p0_gev, 1e-4), 1e-4));

  if (pAxis.size() < 2u || pAxis.size() != logL.size() || !(p0_gev > 0.) || !std::isfinite(p0_gev) || !(span > 0.))
    return fallbackW;

  const double pLo = pAxis.front();
  const double pHi = pAxis.back();
  // Step ~8% of p0, bounded so p0±h stay in the interior away from clamped plateaus.
  double h = 0.08 * std::max(p0_gev, 1e-4);
  h = std::max(h, medianStep * 2.0);
  h = std::min(h, 0.25 * span);
  h = std::min(h, 0.45 * (p0_gev - pLo));
  h = std::min(h, 0.45 * (pHi - p0_gev));
  if (!(h > 0.) || h < medianStep * 0.5) {
    return SigmaFromDeltaNLLHalf(pAxis, logL, p0_gev, fallbackW);
  }

  const double Nm = NLLFromLogLikelihoodCurve(pAxis, logL, p0_gev);
  const double Np = NLLFromLogLikelihoodCurve(pAxis, logL, p0_gev + h);
  const double Nn = NLLFromLogLikelihoodCurve(pAxis, logL, p0_gev - h);
  if (!std::isfinite(Nm) || !std::isfinite(Np) || !std::isfinite(Nn)) return SigmaFromDeltaNLLHalf(pAxis, logL, p0_gev, fallbackW);

  const double d2 = (Np - 2.0 * Nm + Nn) / (h * h);
  if (d2 > 1e-12 && std::isfinite(d2)) {
    const double sig = 1.0 / std::sqrt(d2);
    if (std::isfinite(sig) && sig > 0.) return std::max(sig, 1e-6);
  }
  return SigmaFromDeltaNLLHalf(pAxis, logL, p0_gev, fallbackW);
}

bool ComputeSigmaMEventGeV(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                           const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                           const TVector3& dir1, const TVector3& dir2, double fallback_sigma_m_gev,
                           JointK0sSigmaMEventDiagnostics& out) {
  constexpr double kPionMassGeV = 0.13957;
  out = JointK0sSigmaMEventDiagnostics();

  if (p1Axis.size() < 2u || logL1.size() != p1Axis.size() || p2Axis.size() < 2u || logL2.size() != p2Axis.size())
    return false;
  if (dir1.Mag2() <= 0. || dir2.Mag2() <= 0.) return false;

  auto argmaxP = [](const std::vector<double>& p, const std::vector<double>& l) -> double {
    size_t im = 0;
    for (size_t i = 1; i < l.size(); ++i) {
      if (l[i] > l[im]) im = i;
    }
    return p[im];
  };

  const double p10 = argmaxP(p1Axis, logL1);
  const double p20 = argmaxP(p2Axis, logL2);

  const TVector3 u1 = dir1.Unit();
  const TVector3 u2 = dir2.Unit();
  const double ct = std::max(-1.0, std::min(1.0, u1.Dot(u2)));
  out.cos_theta = ct;

  double m = 0.;
  double d1 = 0.;
  double d2 = 0.;
  if (!PionPiPiInvariantMassAndDerivativesGeV(p10, p20, ct, kPionMassGeV, m, d1, d2)) {
    out.sigma_m_event_gev = fallback_sigma_m_gev;
    return true;
  }
  out.dm_dp1 = d1;
  out.dm_dp2 = d2;

  double sp1 = EstimateSigmaPFromLogLikelihoodCurve(p1Axis, logL1, p10);
  double sp2 = EstimateSigmaPFromLogLikelihoodCurve(p2Axis, logL2, p20);
  if (!(sp1 > 0.) || !std::isfinite(sp1)) sp1 = 0.05 * std::max(p10, 1e-4);
  if (!(sp2 > 0.) || !std::isfinite(sp2)) sp2 = 0.05 * std::max(p20, 1e-4);
  out.sigma_p1_gev = sp1;
  out.sigma_p2_gev = sp2;

  const double sm2 = (d1 * sp1) * (d1 * sp1) + (d2 * sp2) * (d2 * sp2);
  if (!std::isfinite(sm2) || sm2 <= 0.) {
    out.sigma_m_event_gev = fallback_sigma_m_gev;
    return true;
  }

  double sm = std::sqrt(sm2);
  if (ct > kCollinearCosCut) sm *= kCollinearInflate;
  sm = std::max(sm, kFloorSigmaMGeV);
  sm = std::min(sm, kCapSigmaMGeV);
  out.sigma_m_event_gev = sm;
  return true;
}

JointK0sPionMomentumGridResult FitJointMomentaOnGrid(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                                                     const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                                                     const TVector3& dir1, const TVector3& dir2, double pMinGeV,
                                                     double pMaxGeV, double pStepGeV, double mK0sGeV, double sigmaMassGeV,
                                                     double penaltyScale, int refineFactor) {
  JointK0sPionMomentumGridResult out;
  if (!(pStepGeV > 0.) || !(pMaxGeV > pMinGeV) || !std::isfinite(pMinGeV) || !std::isfinite(pMaxGeV)) return out;

  double bestS = -std::numeric_limits<double>::infinity();
  double bestP1 = std::numeric_limits<double>::quiet_NaN();
  double bestP2 = std::numeric_limits<double>::quiet_NaN();
  double bestM = std::numeric_limits<double>::quiet_NaN();

  ScanRectangle(p1Axis, logL1, p2Axis, logL2, dir1, dir2, pMinGeV, pMaxGeV, pMinGeV, pMaxGeV, pStepGeV, mK0sGeV,
                sigmaMassGeV, penaltyScale, bestS, bestP1, bestP2, bestM);

  if (refineFactor >= 2 && std::isfinite(bestS) && std::isfinite(bestP1) && std::isfinite(bestP2)) {
    const double span = 2.5 * pStepGeV;
    const double fineStep = pStepGeV / static_cast<double>(refineFactor);
    const double lo1 = std::max(pMinGeV, bestP1 - span);
    const double hi1 = std::min(pMaxGeV, bestP1 + span);
    const double lo2 = std::max(pMinGeV, bestP2 - span);
    const double hi2 = std::min(pMaxGeV, bestP2 + span);
    ScanRectangle(p1Axis, logL1, p2Axis, logL2, dir1, dir2, lo1, hi1, lo2, hi2, fineStep, mK0sGeV, sigmaMassGeV,
                  penaltyScale, bestS, bestP1, bestP2, bestM);
  }

  if (!std::isfinite(bestS) || bestS <= -0.5 * std::numeric_limits<double>::max()) return out;

  // Post-fit diagnostics only (do not feed back into grid / best point).
  const double logL1b = InterpolateLogLikelihoodClamped(p1Axis, logL1, bestP1);
  const double logL2b = InterpolateLogLikelihoodClamped(p2Axis, logL2, bestP2);
  const double pen = MassPenaltyChi2(bestM, mK0sGeV, sigmaMassGeV, penaltyScale);
  const double c_mass =
      (penaltyScale > 0. && std::isfinite(penaltyScale)) ? (-pen / penaltyScale) : 0.;  // = −(M−m_K0)²/(2σ²)
  const double l_tracks = logL1b + logL2b;
  constexpr double kEpsL = 1e-8;
  constexpr float kRSentinel = 1e6f;
  float rDiag = kRSentinel;
  if (std::isfinite(l_tracks) && std::abs(l_tracks) >= kEpsL && std::isfinite(c_mass)) {
    rDiag = static_cast<float>(std::abs(c_mass) / std::max(std::abs(l_tracks), kEpsL));
  } else if (std::isfinite(l_tracks) && std::abs(l_tracks) < kEpsL) {
    rDiag = kRSentinel;
  } else if (std::isfinite(c_mass)) {
    rDiag = kRSentinel;
  } else {
    rDiag = -999.f;
  }

  out.ok = true;
  out.p1 = static_cast<Float_t>(bestP1);
  out.p2 = static_cast<Float_t>(bestP2);
  out.bestScore = static_cast<Float_t>(bestS);
  out.invMassAtBest = static_cast<Float_t>(bestM);
  out.constraintRatioR = rDiag;
  return out;
}

} // namespace pdJointK0sPionMomentum
