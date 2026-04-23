#include "pdJointK0sPionMomentum.hxx"
#include <TH2F.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace pdJointK0sPionMomentum {
namespace {

constexpr double kCollinearCosCut = 0.95;
constexpr double kCollinearInflate = 2.0;
constexpr double kPionMassMeV = 139.57;
constexpr double kHighlandConstantMeV = 13.6;

double ClampCosine(double c) {
  return std::max(-1.0, std::min(1.0, c));
}

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

bool BuildPionMcsScatteringSamples(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                  std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                  std::vector<double>* rrMidCm) {
  deltaTheta.clear();
  xOverX0.clear();
  if (rrMidCm) rrMidCm->clear();

  std::vector<TVector3> orderedPoints;
  std::vector<double> orderedRR;
  {
    struct HitPoint {
      double rr = 0.0;
      TVector3 pos;
    };
    std::vector<HitPoint> hitPoints;
    int bestPlane = -1;
    size_t bestCount = 0u;
    for (int ipl = 0; ipl < 3; ++ipl) {
      size_t nValid = 0u;
      for (const auto& h : track.Hits[ipl]) {
        if (!std::isfinite(h.Position.X()) || !std::isfinite(h.Position.Y()) || !std::isfinite(h.Position.Z())) continue;
        if (!std::isfinite(static_cast<double>(h.ResidualRange)) || h.ResidualRange < -900.f) continue;
        ++nValid;
      }
      if (nValid > bestCount) {
        bestCount = nValid;
        bestPlane = ipl;
      }
    }
    if (bestPlane >= 0 && bestCount >= 3u) {
      hitPoints.reserve(bestCount);
      for (const auto& h : track.Hits[bestPlane]) {
        if (!std::isfinite(h.Position.X()) || !std::isfinite(h.Position.Y()) || !std::isfinite(h.Position.Z())) continue;
        if (!std::isfinite(static_cast<double>(h.ResidualRange)) || h.ResidualRange < -900.f) continue;
        HitPoint hp;
        hp.rr = static_cast<double>(h.ResidualRange);
        hp.pos = h.Position;
        hitPoints.push_back(hp);
      }
      std::stable_sort(hitPoints.begin(), hitPoints.end(),
                       [](const HitPoint& a, const HitPoint& b) { return a.rr > b.rr; });
      orderedPoints.reserve(hitPoints.size());
      orderedRR.reserve(hitPoints.size());
      for (const auto& hp : hitPoints) {
        orderedPoints.push_back(hp.pos);
        orderedRR.push_back(hp.rr);
      }
    }
  }

  if (orderedPoints.size() < 3u) {
    orderedPoints.clear();
    orderedRR.clear();
    orderedPoints.reserve(track.TrjPoints.size());
    for (const auto& tp : track.TrjPoints) {
      if (!std::isfinite(tp.Position.X()) || !std::isfinite(tp.Position.Y()) || !std::isfinite(tp.Position.Z())) continue;
      orderedPoints.push_back(tp.Position);
    }
    if (orderedPoints.size() >= 2u) {
      std::vector<double> s(orderedPoints.size(), 0.0);
      for (size_t i = 1; i < orderedPoints.size(); ++i) {
        const double dl = (orderedPoints[i] - orderedPoints[i - 1]).Mag();
        s[i] = s[i - 1] + ((std::isfinite(dl) && dl > 0.0) ? dl : 0.0);
      }
      const double totalLen = s.back();
      orderedRR.resize(orderedPoints.size(), 0.0);
      for (size_t i = 0; i < orderedPoints.size(); ++i) orderedRR[i] = std::max(0.0, totalLen - s[i]);
    }
  }
  if (orderedPoints.size() < 3u || orderedRR.size() != orderedPoints.size()) return false;

  const double x0 = (std::isfinite(cfg.radiationLengthCm) && cfg.radiationLengthCm > 1e-9) ? cfg.radiationLengthCm : 14.0;
  const double minSeg = (std::isfinite(cfg.minSegmentLengthCm) && cfg.minSegmentLengthCm > 0.0) ? cfg.minSegmentLengthCm : 0.5;
  const double maxAbsTheta =
      (std::isfinite(cfg.maxAbsDeltaThetaRad) && cfg.maxAbsDeltaThetaRad > 0.0) ? cfg.maxAbsDeltaThetaRad : -1.0;

  for (size_t i = 0; i + 2 < orderedPoints.size(); ++i) {
    const TVector3 d1 = orderedPoints[i + 1] - orderedPoints[i];
    const TVector3 d2 = orderedPoints[i + 2] - orderedPoints[i + 1];
    const double l1 = d1.Mag();
    const double l2 = d2.Mag();
    if (!std::isfinite(l1) || !std::isfinite(l2) || l1 <= 1e-9 || l2 <= 1e-9) continue;
    const double seg = 0.5 * (l1 + l2);
    if (!std::isfinite(seg) || seg < minSeg) continue;
    const double c = ClampCosine(d1.Unit().Dot(d2.Unit()));
    double dTh = std::acos(c);
    if (!std::isfinite(dTh)) continue;
    if (maxAbsTheta > 0.0 && dTh > maxAbsTheta) dTh = maxAbsTheta;
    const double xov = seg / x0;
    if (!std::isfinite(xov) || xov <= 0.0) continue;
    const double rr = orderedRR[i + 1];
    if (!std::isfinite(rr) || rr < 0.0) continue;
    deltaTheta.push_back(dTh);
    xOverX0.push_back(xov);
    if (rrMidCm) rrMidCm->push_back(rr);
  }
  if (deltaTheta.empty() || deltaTheta.size() != xOverX0.size()) return false;
  if (rrMidCm && rrMidCm->size() != deltaTheta.size()) return false;
  return true;
}

MCSLikelihood::MCSLikelihood(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg) {
  theta0_floor_rad_ = (std::isfinite(cfg.theta0FloorRad) && cfg.theta0FloorRad > 0.) ? cfg.theta0FloorRad : 1e-6;
  BuildPionMcsScatteringSamples(track, cfg, delta_theta_, x_over_x0_, nullptr);
}

double MCSLikelihood::ComputeNLL(double momentumGeV) const {
  if (delta_theta_.empty() || x_over_x0_.size() != delta_theta_.size()) return std::numeric_limits<double>::infinity();
  if (!std::isfinite(momentumGeV) || momentumGeV <= 0.) return std::numeric_limits<double>::infinity();

  const double pMeV = momentumGeV * 1000.0;
  const double eMeV = std::sqrt(pMeV * pMeV + kPionMassMeV * kPionMassMeV);
  if (!std::isfinite(eMeV) || eMeV <= 0.) return std::numeric_limits<double>::infinity();
  const double beta = pMeV / eMeV;
  if (!std::isfinite(beta) || beta <= 0.) return std::numeric_limits<double>::infinity();

  double nll = 0.;
  for (size_t i = 0; i < delta_theta_.size(); ++i) {
    const double xOverX0 = x_over_x0_[i];
    const double logArg = std::max(xOverX0, 1e-12);
    double corr = 1.0 + 0.038 * std::log(logArg);
    if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;

    double theta0 = (kHighlandConstantMeV / (beta * pMeV)) * std::sqrt(xOverX0) * corr;
    if (!std::isfinite(theta0) || theta0 <= 0.) theta0 = theta0_floor_rad_;
    theta0 = std::max(theta0, theta0_floor_rad_);

    const double dt = delta_theta_[i];
    const double pull = dt / theta0;
    nll += 0.5 * pull * pull + std::log(theta0);
  }
  return (std::isfinite(nll)) ? nll : std::numeric_limits<double>::infinity();
}

bool BuildPionMCSLogLikelihoodVsMomentumCurve(const AnaParticlePD& track, const std::vector<double>& pAxisGeV,
                                              const MCSLikelihoodConfig& cfg, std::vector<double>& logL) {
  logL.clear();
  if (pAxisGeV.empty()) return false;
  const MCSLikelihood mcs(track, cfg);
  if (!mcs.HasSamples()) return false;
  logL.reserve(pAxisGeV.size());
  for (size_t i = 0; i < pAxisGeV.size(); ++i) {
    const double nll = mcs.ComputeNLL(pAxisGeV[i]);
    if (!std::isfinite(nll)) return false;
    logL.push_back(-nll);
  }
  return logL.size() == pAxisGeV.size();
}

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
                           JointK0sSigmaMEventDiagnostics& out, double sigmaMinGeV, double sigmaMaxGeV) {
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
  const double sMin = std::isfinite(sigmaMinGeV) ? sigmaMinGeV : 0.005;
  const double sMax = std::isfinite(sigmaMaxGeV) ? sigmaMaxGeV : 0.05;
  const double lo = std::min(sMin, sMax);
  const double hi = std::max(sMin, sMax);
  sm = std::max(sm, lo);
  sm = std::min(sm, hi);
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

bool MakeJointK0sObjectiveTH2CoarsePass(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                                        const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                                        const TVector3& dir1, const TVector3& dir2, double pMinGeV, double pMaxGeV,
                                        double pStepGeV, double mK0sGeV, double sigmaMassGeV, double penaltyScale,
                                        const char* nameObjective, const char* titleObjective, const char* namePenalty,
                                        const char* titlePenalty, const char* nameTrackLL, const char* titleTrackLL,
                                        TH2F*& hObjective, TH2F*& hMassPenalty, TH2F*& hTrackLogLSum) {
  hObjective = nullptr;
  hMassPenalty = nullptr;
  hTrackLogLSum = nullptr;
  if (!nameObjective || !namePenalty || !nameTrackLL || !titleObjective || !titlePenalty || !titleTrackLL)
    return false;
  if (!(pStepGeV > 0.) || !(pMaxGeV > pMinGeV) || !std::isfinite(pMinGeV) || !std::isfinite(pMaxGeV)) return false;
  const int nSteps = static_cast<int>(std::floor((pMaxGeV - pMinGeV) / pStepGeV + 1e-9)) + 1;
  if (nSteps < 2) return false;
  const double xLo = pMinGeV - 0.5 * pStepGeV;
  const double xHi = pMaxGeV + 0.5 * pStepGeV;
  TH2F* hS = new TH2F(nameObjective, titleObjective, nSteps, xLo, xHi, nSteps, xLo, xHi);
  TH2F* hP = new TH2F(namePenalty, titlePenalty, nSteps, xLo, xHi, nSteps, xLo, xHi);
  TH2F* hLL = new TH2F(nameTrackLL, titleTrackLL, nSteps, xLo, xHi, nSteps, xLo, xHi);
  if (!hS || !hP || !hLL) {
    delete hS;
    delete hP;
    delete hLL;
    return false;
  }
  hS->SetStats(0);
  hP->SetStats(0);
  hLL->SetStats(0);

  auto massPenaltyChi2 = [&](double massGeV) -> double {
    if (!(penaltyScale > 0.) || !std::isfinite(penaltyScale)) return 0.;
    const double sig = (sigmaMassGeV > 1e-9 && std::isfinite(sigmaMassGeV)) ? sigmaMassGeV : 1e-3;
    if (!std::isfinite(massGeV) || !std::isfinite(mK0sGeV)) return 1.e9;
    const double d = (massGeV - mK0sGeV) / sig;
    return 0.5 * penaltyScale * d * d;
  };

  int nFilled = 0;
  for (double a1 = pMinGeV; a1 <= pMaxGeV + 1e-12; a1 += pStepGeV) {
    const double L1 = InterpolateLogLikelihoodClamped(p1Axis, logL1, a1);
    if (!std::isfinite(L1)) continue;
    for (double a2 = pMinGeV; a2 <= pMaxGeV + 1e-12; a2 += pStepGeV) {
      const double L2 = InterpolateLogLikelihoodClamped(p2Axis, logL2, a2);
      if (!std::isfinite(L2)) continue;
      double mPiPi = 0.;
      if (!PionPairInvariantMassGeV(a1, a2, dir1, dir2, mPiPi)) continue;
      const double pen = massPenaltyChi2(mPiPi);
      const double S = L1 + L2 - pen;
      if (!std::isfinite(S) || !std::isfinite(pen)) continue;
      const int ix = hS->GetXaxis()->FindBin(a1);
      const int iy = hS->GetYaxis()->FindBin(a2);
      if (ix < 1 || iy < 1 || ix > hS->GetNbinsX() || iy > hS->GetNbinsY()) continue;
      hS->SetBinContent(ix, iy, S);
      hP->SetBinContent(ix, iy, pen);
      hLL->SetBinContent(ix, iy, L1 + L2);
      ++nFilled;
    }
  }
  if (nFilled == 0) {
    delete hS;
    delete hP;
    delete hLL;
    return false;
  }
  hObjective = hS;
  hMassPenalty = hP;
  hTrackLogLSum = hLL;
  return true;
}

} // namespace pdJointK0sPionMomentum
