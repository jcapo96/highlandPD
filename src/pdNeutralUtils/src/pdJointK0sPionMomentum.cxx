#include "pdJointK0sPionMomentum.hxx"
#include <TH2F.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace pdJointK0sPionMomentum {
namespace {

constexpr double kCollinearCosCut = 0.95;
constexpr double kCollinearInflate = 2.0;
constexpr double kHighlandConstantMeV = 13.6;

inline float BetheBloch_dEdx_MeV_per_cm(float beta, float gamma) {
  // Constants for liquid argon
  const float K = 0.307075f; // MeV mol^-1 cm^2
  const float Z = 18.0f;
  const float A = 39.948f;   // g/mol
  const float rho = 1.396f;  // g/cm^3
  const float I = 188e-6f;   // MeV

  const float me = 0.511f;   // MeV
  const float z = 1.0f;      // pion charge

  float beta2 = beta * beta;
  if (beta2 < 1e-6f) beta2 = 1e-6f;

  const float gamma2 = gamma * gamma;

  const float Wmax = (2.0f * me * beta2 * gamma2);

  float argument = (2.0f * me * beta2 * gamma2 * Wmax) / (I * I);
  argument = std::max(argument, 1.0001f);

  const float logTerm = 0.5f * std::log(argument);

  const float dEdx = K * (Z / A) * (z * z / beta2) * (logTerm - beta2);

  return dEdx * rho; // MeV/cm
}

bool HasValidXYZ(const TVector3& p) {
  return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900.0 && p.Y() > -900.0 &&
         p.Z() > -900.0;
}

bool BuildStableLocalFrame(const TVector3& zAxis, TVector3& xAxis, TVector3& yAxis) {
  if (!(zAxis.Mag2() > 1e-12)) return false;
  const TVector3 z = zAxis.Unit();
  TVector3 ref(1.0, 0.0, 0.0);
  if (std::abs(z.Dot(ref)) > 0.9) ref.SetXYZ(0.0, 1.0, 0.0);
  xAxis = ref.Cross(z);
  if (!(xAxis.Mag2() > 1e-12)) return false;
  xAxis = xAxis.Unit();
  yAxis = z.Cross(xAxis);
  if (!(yAxis.Mag2() > 1e-12)) return false;
  yAxis = yAxis.Unit();
  return true;
}

bool FitDirectionFromPoints(const std::vector<TVector3>& pts, TVector3& dirOut, double* chi2NdfOut = nullptr) {
  if (pts.size() < 2u) return false;
  double mx = 0.0, my = 0.0, mz = 0.0;
  for (const auto& p : pts) {
    mx += p.X();
    my += p.Y();
    mz += p.Z();
  }
  const double invN = 1.0 / static_cast<double>(pts.size());
  mx *= invN;
  my *= invN;
  mz *= invN;

  double cxx = 0.0, cxy = 0.0, cxz = 0.0, cyy = 0.0, cyz = 0.0, czz = 0.0;
  for (const auto& p : pts) {
    const double dx = p.X() - mx;
    const double dy = p.Y() - my;
    const double dz = p.Z() - mz;
    cxx += dx * dx;
    cxy += dx * dy;
    cxz += dx * dz;
    cyy += dy * dy;
    cyz += dy * dz;
    czz += dz * dz;
  }
  TVector3 v(0.0, 0.0, 1.0);
  for (int it = 0; it < 12; ++it) {
    TVector3 w(cxx * v.X() + cxy * v.Y() + cxz * v.Z(), cxy * v.X() + cyy * v.Y() + cyz * v.Z(),
               cxz * v.X() + cyz * v.Y() + czz * v.Z());
    if (!(w.Mag2() > 1e-20)) break;
    v = w.Unit();
  }
  if (!(v.Mag2() > 1e-12)) return false;
  dirOut = v.Unit();
  if (chi2NdfOut) {
    const TVector3 c(mx, my, mz);
    double chi2 = 0.0;
    for (const auto& p : pts) {
      const TVector3 d = p - c;
      const double par = d.Dot(dirOut);
      const TVector3 perp = d - par * dirOut;
      chi2 += perp.Mag2();
    }
    const int ndf = std::max(1, static_cast<int>(pts.size()) - 2);
    *chi2NdfOut = chi2 / static_cast<double>(ndf);
  }
  return true;
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

struct McsInternalSegment {
  TVector3 start;
  TVector3 end;
  TVector3 direction;
  std::vector<TVector3> points;
  double arcLengthCm = 0.0;
  double rrMidCm = 0.0;
  double fitChi2Ndf = -999.0;
};

/// Shared segmentation: build the per-segment list used by the MCS estimator.
/// Returns false if input is too short or geometry is degenerate.
bool BuildMcsInternalSegments(const std::vector<TVector3>& orderedPoints, const MCSLikelihoodConfig& cfg,
                              std::vector<McsInternalSegment>& segments) {
  segments.clear();
  if (orderedPoints.size() < 2u) return false;

  std::vector<double> arcLen(orderedPoints.size(), 0.0);
  for (size_t i = 1; i < orderedPoints.size(); ++i) {
    const double dl = (orderedPoints[i] - orderedPoints[i - 1]).Mag();
    if (!std::isfinite(dl) || dl <= 0.0) continue;
    arcLen[i] = arcLen[i - 1] + dl;
  }
  const double totalLen = arcLen.back();
  if (!std::isfinite(totalLen) || totalLen <= 0.0) return false;

  std::vector<double> orderedRR(orderedPoints.size(), 0.0);
  for (size_t i = 0; i < orderedPoints.size(); ++i) orderedRR[i] = std::max(0.0, totalLen - arcLen[i]);

  const double minSeg = (std::isfinite(cfg.minSegmentLengthCm) && cfg.minSegmentLengthCm > 0.0) ? cfg.minSegmentLengthCm : 0.5;
  const double targetSeg =
      (std::isfinite(cfg.targetSegmentLengthCm) && cfg.targetSegmentLengthCm > 0.0) ? cfg.targetSegmentLengthCm : 10.0;

  segments.reserve(orderedPoints.size() / 2u);
  size_t segStart = 0u;
  double accum = 0.0;
  for (size_t i = 1; i < orderedPoints.size(); ++i) {
    const double dl = (orderedPoints[i] - orderedPoints[i - 1]).Mag();
    if (std::isfinite(dl) && dl > 0.0) accum += dl;
    if (accum < targetSeg && i + 1 < orderedPoints.size()) continue;

    const TVector3 disp = orderedPoints[i] - orderedPoints[segStart];
    const double chord = disp.Mag();
    if (std::isfinite(chord) && chord > 1e-9 && std::isfinite(accum) && accum >= minSeg) {
      McsInternalSegment seg;
      seg.start = orderedPoints[segStart];
      seg.end = orderedPoints[i];
      seg.points.reserve(i - segStart + 1u);
      for (size_t j = segStart; j <= i; ++j) seg.points.push_back(orderedPoints[j]);
      TVector3 fitDir;
      double fitChi2Ndf = -999.0;
      if (FitDirectionFromPoints(seg.points, fitDir, &fitChi2Ndf)) {
        if (fitDir.Dot(disp) < 0.0) fitDir *= -1.0;
        seg.direction = fitDir;
        seg.fitChi2Ndf = fitChi2Ndf;
      } else {
        seg.direction = disp.Unit();
      }
      seg.arcLengthCm = accum;
      seg.rrMidCm = 0.5 * (orderedRR[segStart] + orderedRR[i]);
      segments.push_back(seg);
    }
    segStart = i;
    accum = 0.0;
  }
  return !segments.empty();
}

bool BuildPionMcsScatteringSamplesCore(const std::vector<TVector3>& orderedPoints, const MCSLikelihoodConfig& cfg,
                                      std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                      std::vector<double>* rrMidCm, std::vector<double>* thetaXZ,
                                      std::vector<double>* thetaYZ, std::vector<double>* segmentFitChi2Ndf,
                                      std::vector<double>* segment1RrToEndCm, std::vector<double>* segment2RrToEndCm,
                                      std::vector<double>* segmentLengthCm, std::vector<double>* segmentFitChi2NdfSingle,
                                      std::vector<double>* segmentRrToEndCmSingle) {
  std::vector<McsInternalSegment> segments;
  if (!BuildMcsInternalSegments(orderedPoints, cfg, segments)) return false;

  const double x0 = (std::isfinite(cfg.radiationLengthCm) && cfg.radiationLengthCm > 1e-9) ? cfg.radiationLengthCm : 14.0;
  const double minSeg = (std::isfinite(cfg.minSegmentLengthCm) && cfg.minSegmentLengthCm > 0.0) ? cfg.minSegmentLengthCm : 0.5;
  const double maxAbsTheta =
      (std::isfinite(cfg.maxAbsDeltaThetaRad) && cfg.maxAbsDeltaThetaRad > 0.0) ? cfg.maxAbsDeltaThetaRad : -1.0;

  if (segments.size() < 2u) return false;
  if (segmentLengthCm || segmentFitChi2NdfSingle || segmentRrToEndCmSingle) {
    for (size_t i = 0; i < segments.size(); ++i) {
      if (segmentLengthCm) segmentLengthCm->push_back(segments[i].arcLengthCm);
      if (segmentFitChi2NdfSingle) segmentFitChi2NdfSingle->push_back(segments[i].fitChi2Ndf);
      if (segmentRrToEndCmSingle) segmentRrToEndCmSingle->push_back(segments[i].rrMidCm);
    }
  }

  for (size_t i = 0; i + 1 < segments.size(); ++i) {
    const double seg = 0.5 * (segments[i].arcLengthCm + segments[i + 1].arcLengthCm);
    if (!std::isfinite(seg) || seg < minSeg) continue;
    TVector3 ex, ey;
    if (!BuildStableLocalFrame(segments[i].direction, ex, ey)) continue;
    const TVector3 z = segments[i].direction.Unit();
    const TVector3 v = segments[i + 1].direction.Unit();
    const double vz = v.Dot(z);
    if (!std::isfinite(vz) || std::abs(vz) < 1e-9) continue;
    const double vx = v.Dot(ex);
    const double vy = v.Dot(ey);
    double dThXZ = std::atan2(vx, vz);
    double dThYZ = std::atan2(vy, vz);
    if (!std::isfinite(dThXZ) || !std::isfinite(dThYZ)) continue;
    if (maxAbsTheta > 0.0) {
      dThXZ = std::max(-maxAbsTheta, std::min(maxAbsTheta, dThXZ));
      dThYZ = std::max(-maxAbsTheta, std::min(maxAbsTheta, dThYZ));
    }
    const double dTh = std::sqrt(dThXZ * dThXZ + dThYZ * dThYZ);
    const double xov = seg / x0;
    if (!std::isfinite(xov) || xov <= 0.0) continue;
    const double rr = 0.5 * (segments[i].rrMidCm + segments[i + 1].rrMidCm);
    if (!std::isfinite(rr) || rr < 0.0) continue;
    deltaTheta.push_back(dTh);
    xOverX0.push_back(xov);
    if (rrMidCm) rrMidCm->push_back(rr);
    if (thetaXZ) thetaXZ->push_back(dThXZ);
    if (thetaYZ) thetaYZ->push_back(dThYZ);
    if (segmentFitChi2Ndf) {
      const double q1 = segments[i].fitChi2Ndf;
      const double q2 = segments[i + 1].fitChi2Ndf;
      if (std::isfinite(q1) && q1 >= 0.0 && std::isfinite(q2) && q2 >= 0.0)
        segmentFitChi2Ndf->push_back(0.5 * (q1 + q2));
      else if (std::isfinite(q1) && q1 >= 0.0)
        segmentFitChi2Ndf->push_back(q1);
      else if (std::isfinite(q2) && q2 >= 0.0)
        segmentFitChi2Ndf->push_back(q2);
      else
        segmentFitChi2Ndf->push_back(-999.0);
    }
    if (segment1RrToEndCm) segment1RrToEndCm->push_back(segments[i].rrMidCm);
    if (segment2RrToEndCm) segment2RrToEndCm->push_back(segments[i + 1].rrMidCm);
  }

  if (deltaTheta.empty() || deltaTheta.size() != xOverX0.size()) return false;
  if (rrMidCm && rrMidCm->size() != deltaTheta.size()) return false;
  if (thetaXZ && thetaXZ->size() != deltaTheta.size()) return false;
  if (thetaYZ && thetaYZ->size() != deltaTheta.size()) return false;
  if (segmentFitChi2Ndf && segmentFitChi2Ndf->size() != deltaTheta.size()) return false;
  if (segment1RrToEndCm && segment1RrToEndCm->size() != deltaTheta.size()) return false;
  if (segment2RrToEndCm && segment2RrToEndCm->size() != deltaTheta.size()) return false;
  if (segmentLengthCm && segmentLengthCm->size() != segments.size()) return false;
  if (segmentFitChi2NdfSingle && segmentFitChi2NdfSingle->size() != segments.size()) return false;
  if (segmentRrToEndCmSingle && segmentRrToEndCmSingle->size() != segments.size()) return false;
  return true;
}

} // namespace

bool BuildPionMcsScatteringSamples(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                  std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                  std::vector<double>* rrMidCm, std::vector<double>* thetaXZ,
                                  std::vector<double>* thetaYZ, std::vector<double>* segmentFitChi2Ndf,
                                  std::vector<double>* segment1RrToEndCm, std::vector<double>* segment2RrToEndCm,
                                  std::vector<double>* segmentLengthCm, std::vector<double>* segmentFitChi2NdfSingle,
                                  std::vector<double>* segmentRrToEndCmSingle) {
  deltaTheta.clear();
  xOverX0.clear();
  if (rrMidCm) rrMidCm->clear();
  if (thetaXZ) thetaXZ->clear();
  if (thetaYZ) thetaYZ->clear();
  if (segmentFitChi2Ndf) segmentFitChi2Ndf->clear();
  if (segment1RrToEndCm) segment1RrToEndCm->clear();
  if (segment2RrToEndCm) segment2RrToEndCm->clear();
  if (segmentLengthCm) segmentLengthCm->clear();
  if (segmentFitChi2NdfSingle) segmentFitChi2NdfSingle->clear();
  if (segmentRrToEndCmSingle) segmentRrToEndCmSingle->clear();

  std::vector<TVector3> orderedPoints;
  orderedPoints.reserve(track.TrjPoints.size());
  for (size_t i = 0; i < track.TrjPoints.size(); ++i) {
    if (i == 0) continue; // first calo sample often off interior trend; still drawn in ED, not used for MCS
    const AnaTrajectoryPointPD& tp = track.TrjPoints[i];
    if (!tp.IsInTPC) continue;
    if (!HasValidXYZ(tp.Position)) continue;
    orderedPoints.push_back(tp.Position);
  }
  if (orderedPoints.size() < 3u) return false;
  return BuildPionMcsScatteringSamplesCore(orderedPoints, cfg, deltaTheta, xOverX0, rrMidCm, thetaXZ, thetaYZ,
                                          segmentFitChi2Ndf, segment1RrToEndCm, segment2RrToEndCm, segmentLengthCm,
                                          segmentFitChi2NdfSingle, segmentRrToEndCmSingle);
}

bool BuildPionMcsScatteringSamplesFromOrderedPositions(const std::vector<TVector3>& orderedPositions,
                                                       const MCSLikelihoodConfig& cfg,
                                                       std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                                       std::vector<double>* rrMidCm, std::vector<double>* thetaXZ,
                                                       std::vector<double>* thetaYZ, std::vector<double>* segmentFitChi2Ndf,
                                                       std::vector<double>* segment1RrToEndCm,
                                                       std::vector<double>* segment2RrToEndCm,
                                                       std::vector<double>* segmentLengthCm,
                                                       std::vector<double>* segmentFitChi2NdfSingle,
                                                       std::vector<double>* segmentRrToEndCmSingle) {
  deltaTheta.clear();
  xOverX0.clear();
  if (rrMidCm) rrMidCm->clear();
  if (thetaXZ) thetaXZ->clear();
  if (thetaYZ) thetaYZ->clear();
  if (segmentFitChi2Ndf) segmentFitChi2Ndf->clear();
  if (segment1RrToEndCm) segment1RrToEndCm->clear();
  if (segment2RrToEndCm) segment2RrToEndCm->clear();
  if (segmentLengthCm) segmentLengthCm->clear();
  if (segmentFitChi2NdfSingle) segmentFitChi2NdfSingle->clear();
  if (segmentRrToEndCmSingle) segmentRrToEndCmSingle->clear();
  if (orderedPositions.size() < 3u) return false;
  return BuildPionMcsScatteringSamplesCore(orderedPositions, cfg, deltaTheta, xOverX0, rrMidCm, thetaXZ, thetaYZ,
                                          segmentFitChi2Ndf, segment1RrToEndCm, segment2RrToEndCm, segmentLengthCm,
                                          segmentFitChi2NdfSingle, segmentRrToEndCmSingle);
}

namespace {
void CopyInternalSegmentsToGeometry(const std::vector<McsInternalSegment>& segs,
                                    std::vector<MCSSegmentGeometry>& out) {
  out.clear();
  out.reserve(segs.size());
  for (const auto& s : segs) {
    MCSSegmentGeometry g;
    g.startPoint = s.start;
    g.endPoint = s.end;
    TVector3 c(0.0, 0.0, 0.0);
    if (!s.points.empty()) {
      for (const auto& p : s.points) c += p;
      c *= (1.0 / static_cast<double>(s.points.size()));
    } else {
      c = 0.5 * (s.start + s.end);
    }
    g.centroid = c;
    g.fittedDirection = (s.direction.Mag2() > 1e-20) ? s.direction.Unit() : (s.end - s.start).Unit();
    g.arcLengthCm = s.arcLengthCm;
    out.push_back(g);
  }
}
} // namespace

bool BuildPionMcsSegmentGeometry(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                 std::vector<MCSSegmentGeometry>& outSegments) {
  outSegments.clear();
  std::vector<TVector3> orderedPoints;
  orderedPoints.reserve(track.TrjPoints.size());
  for (size_t i = 0; i < track.TrjPoints.size(); ++i) {
    if (i == 0) continue; // match BuildPionMcsScatteringSamples: skip first TrjPoints entry for MCS
    const AnaTrajectoryPointPD& tp = track.TrjPoints[i];
    if (!tp.IsInTPC) continue;
    if (!HasValidXYZ(tp.Position)) continue;
    orderedPoints.push_back(tp.Position);
  }
  if (orderedPoints.size() < 3u) return false;
  std::vector<McsInternalSegment> segs;
  if (!BuildMcsInternalSegments(orderedPoints, cfg, segs)) return false;
  CopyInternalSegmentsToGeometry(segs, outSegments);
  return !outSegments.empty();
}

bool BuildPionMcsSegmentGeometryFromOrderedPositions(const std::vector<TVector3>& orderedPositions,
                                                      const MCSLikelihoodConfig& cfg,
                                                      std::vector<MCSSegmentGeometry>& outSegments) {
  outSegments.clear();
  if (orderedPositions.size() < 3u) return false;
  std::vector<McsInternalSegment> segs;
  if (!BuildMcsInternalSegments(orderedPositions, cfg, segs)) return false;
  CopyInternalSegmentsToGeometry(segs, outSegments);
  return !outSegments.empty();
}

MCSLikelihood::MCSLikelihood(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg) {
  theta0_floor_rad_ = (std::isfinite(cfg.theta0FloorRad) && cfg.theta0FloorRad > 0.) ? cfg.theta0FloorRad : 1e-6;
  radiation_length_cm_ =
      (std::isfinite(cfg.radiationLengthCm) && cfg.radiationLengthCm > 1e-9) ? cfg.radiationLengthCm : 14.0;
  use_detector_sigma_ = cfg.useDetectorSigma;
  detector_sigma_floor_rad_ =
      (std::isfinite(cfg.detectorSigmaFloorRad) && cfg.detectorSigmaFloorRad > 0.) ? cfg.detectorSigmaFloorRad : 1e-6;
  detector_sigma_a_ = (std::isfinite(cfg.detectorSigmaA) && cfg.detectorSigmaA >= 0.) ? cfg.detectorSigmaA : 0.0;
  detector_sigma_c_ = std::isfinite(cfg.detectorSigmaC) ? cfg.detectorSigmaC : 0.0;
  if (use_detector_sigma_) {
    // Keep only analytic detector-resolution model:
    // sigma_det^2(L) = A/L^2 + C, with floor protection.
    if (!(detector_sigma_a_ > 0.) && !(detector_sigma_c_ > 0.)) use_detector_sigma_ = false;
  }
  std::vector<double> deltaThetaUnused;
  BuildPionMcsScatteringSamples(track, cfg, deltaThetaUnused, x_over_x0_, nullptr, &theta_xz_, &theta_yz_);
  cumulative_length_cm_.clear();
  cumulative_length_cm_.reserve(x_over_x0_.size());
  double cumulativeLength = 0.0;
  for (size_t i = 0; i < x_over_x0_.size(); ++i) {
    const double segmentLengthCm = x_over_x0_[i] * radiation_length_cm_;
    const double midpoint = cumulativeLength + 0.5 * segmentLengthCm;
    cumulative_length_cm_.push_back(midpoint);
    cumulativeLength += segmentLengthCm;
  }
}

double MCSLikelihood::ComputeNLL(double momentumGeV) const {
  if (theta_xz_.empty() || theta_yz_.size() != theta_xz_.size() || x_over_x0_.size() != theta_xz_.size() ||
      cumulative_length_cm_.size() != theta_xz_.size())
    return std::numeric_limits<double>::infinity();
  if (!std::isfinite(momentumGeV) || momentumGeV <= 0.) return std::numeric_limits<double>::infinity();

  const float m_pi = 139.57f; // MeV
  const float p0_MeV = static_cast<float>(momentumGeV * 1000.0);
  const float E0 = std::sqrt(p0_MeV * p0_MeV + m_pi * m_pi);
  if (!std::isfinite(E0) || E0 <= 0.0f) return std::numeric_limits<double>::infinity();

  double nll = 0.;
  for (size_t i = 0; i < theta_xz_.size(); ++i) {
    float Ej = E0;
    const float step = 1.0f; // cm
    const float L = static_cast<float>(cumulative_length_cm_[i]);
    for (float x = 0.0f; x < L; x += step) {
      const float p2 = Ej * Ej - m_pi * m_pi;
      if (!std::isfinite(p2) || p2 <= 0.0f) {
        Ej = m_pi;
        break;
      }
      const float pj = std::sqrt(p2);
      const float betaj = pj / Ej;
      const float gammaj = Ej / m_pi;
      const float dEdx = BetheBloch_dEdx_MeV_per_cm(betaj, gammaj);
      Ej -= dEdx * step;
      if (Ej <= m_pi + 1.0f) break;
    }
    if (Ej <= m_pi + 1.0f) continue;

    const float p2local = Ej * Ej - m_pi * m_pi;
    if (!std::isfinite(p2local) || p2local <= 0.0f) continue;
    const float pj = std::sqrt(p2local);
    const float betaj = pj / Ej;
    if (!std::isfinite(betaj) || betaj <= 0.0f || !std::isfinite(pj) || pj <= 0.0f) continue;

    const double xOverX0 = x_over_x0_[i];
    const double logArg = std::max(xOverX0, 1e-12);
    double corr = 1.0 + 0.038 * std::log(logArg);
    if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;

    double theta0 = (kHighlandConstantMeV / (static_cast<double>(betaj) * static_cast<double>(pj))) *
                    std::sqrt(xOverX0) * corr;
    if (!std::isfinite(theta0) || theta0 <= 0.) theta0 = theta0_floor_rad_;
    theta0 = std::max(theta0, theta0_floor_rad_);

    double sigma = theta0;
    if (use_detector_sigma_) {
      const double segCm = xOverX0 * radiation_length_cm_;
      double sigmaDet2 = 0.0;
      if (std::isfinite(segCm) && segCm > 0.0) sigmaDet2 = detector_sigma_a_ / (segCm * segCm) + detector_sigma_c_;
      if (!std::isfinite(sigmaDet2) || sigmaDet2 <= 0.0) sigmaDet2 = detector_sigma_floor_rad_ * detector_sigma_floor_rad_;
      double sigmaDet = std::sqrt(sigmaDet2);
      sigmaDet = std::max(sigmaDet, detector_sigma_floor_rad_);
      sigma = std::sqrt(theta0 * theta0 + sigmaDet * sigmaDet);
    }
    sigma = std::max(sigma, theta0_floor_rad_);

    const double dtxz = theta_xz_[i];
    const double dtyz = theta_yz_[i];
    if (!std::isfinite(dtxz) || !std::isfinite(dtyz)) return std::numeric_limits<double>::infinity();
    const double pullX = dtxz / sigma;
    const double pullY = dtyz / sigma;
    nll += 0.5 * (pullX * pullX + pullY * pullY) + 2.0 * std::log(sigma);
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
