#include "pdMomReconstructionMCS.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

constexpr double kHighlandConstantMeV = 13.6;

inline float BetheBloch_dEdx_MeV_per_cm(float beta, float gamma) {
  const float K = 0.307075f;
  const float Z = 18.0f;
  const float A = 39.948f;
  const float rho = 1.396f;
  const float I = 188e-6f;
  const float me = 0.511f;
  const float z = 1.0f;

  float beta2 = beta * beta;
  if (beta2 < 1e-6f) beta2 = 1e-6f;
  const float gamma2 = gamma * gamma;
  const float Wmax = (2.0f * me * beta2 * gamma2);
  float argument = (2.0f * me * beta2 * gamma2 * Wmax) / (I * I);
  argument = std::max(argument, 1.0001f);
  const float logTerm = 0.5f * std::log(argument);
  const float dEdx = K * (Z / A) * (z * z / beta2) * (logTerm - beta2);
  return dEdx * rho;
}

void FillMomentumScanAxis(double pMcsMaxGeV, std::vector<double>& pAxisGeV) {
  pAxisGeV.clear();
  const double pMcsMax = (std::isfinite(pMcsMaxGeV) && pMcsMaxGeV > 0.05) ? pMcsMaxGeV : 2.5;
  for (double p = 0.05; p <= pMcsMax + 1e-12; p += 0.01) pAxisGeV.push_back(p);
}

}  // namespace

namespace pdMomReconstruction {
namespace {

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

struct McsInternalSegment {
  TVector3 start;
  TVector3 end;
  TVector3 direction;
  std::vector<TVector3> points;
  double arcLengthCm = 0.0;
  double rrMidCm = 0.0;
  double fitChi2Ndf = -999.0;
};

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

void CopyInternalSegmentsToGeometry(const std::vector<McsInternalSegment>& segs, std::vector<MCSSegmentGeometry>& out) {
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

}  // namespace

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
    if (i == 0) continue;
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

bool BuildPionMcsSegmentGeometry(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                 std::vector<MCSSegmentGeometry>& outSegments) {
  outSegments.clear();
  std::vector<TVector3> orderedPoints;
  orderedPoints.reserve(track.TrjPoints.size());
  for (size_t i = 0; i < track.TrjPoints.size(); ++i) {
    if (i == 0) continue;
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

double ComputePionMcsNLLFromScatteringSamples(const std::vector<double>& theta_xz, const std::vector<double>& theta_yz,
                                              const std::vector<double>& x_over_x0,
                                              const std::vector<double>& cumulative_length_cm_midtrack,
                                              double momentumGeV, const MCSLikelihoodConfig& cfg) {
  if (theta_xz.empty() || theta_yz.size() != theta_xz.size() || x_over_x0.size() != theta_xz.size() ||
      cumulative_length_cm_midtrack.size() != theta_xz.size())
    return std::numeric_limits<double>::infinity();
  if (!std::isfinite(momentumGeV) || momentumGeV <= 0.) return std::numeric_limits<double>::infinity();

  double theta0_floor_rad = (std::isfinite(cfg.theta0FloorRad) && cfg.theta0FloorRad > 0.) ? cfg.theta0FloorRad : 1e-6;
  double radiation_length_cm =
      (std::isfinite(cfg.radiationLengthCm) && cfg.radiationLengthCm > 1e-9) ? cfg.radiationLengthCm : 14.0;
  bool use_detector_sigma = cfg.useDetectorSigma;
  double detector_sigma_floor_rad =
      (std::isfinite(cfg.detectorSigmaFloorRad) && cfg.detectorSigmaFloorRad > 0.) ? cfg.detectorSigmaFloorRad : 1e-6;
  double detector_sigma_a = (std::isfinite(cfg.detectorSigmaA) && cfg.detectorSigmaA >= 0.) ? cfg.detectorSigmaA : 0.0;
  double detector_sigma_c = std::isfinite(cfg.detectorSigmaC) ? cfg.detectorSigmaC : 0.0;
  if (use_detector_sigma) {
    if (!(detector_sigma_a > 0.) && !(detector_sigma_c > 0.)) use_detector_sigma = false;
  }

  const float m_pi = 139.57f;
  const float p0_MeV = static_cast<float>(momentumGeV * 1000.0);
  const float E0 = std::sqrt(p0_MeV * p0_MeV + m_pi * m_pi);
  if (!std::isfinite(E0) || E0 <= 0.0f) return std::numeric_limits<double>::infinity();

  double nll = 0.;
  for (size_t i = 0; i < theta_xz.size(); ++i) {
    float Ej = E0;
    const float step = 1.0f;
    const float L = static_cast<float>(cumulative_length_cm_midtrack[i]);
    if (!std::isfinite(L) || L < 0.0f) return std::numeric_limits<double>::infinity();
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

    const double xOverX0 = x_over_x0[i];
    if (!std::isfinite(xOverX0) || xOverX0 <= 0.0) return std::numeric_limits<double>::infinity();
    const double logArg = std::max(xOverX0, 1e-12);
    double corr = 1.0 + 0.038 * std::log(logArg);
    if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;

    double theta0 = (kHighlandConstantMeV / (static_cast<double>(betaj) * static_cast<double>(pj))) * std::sqrt(xOverX0) *
                    corr;
    if (!std::isfinite(theta0) || theta0 <= 0.) theta0 = theta0_floor_rad;
    theta0 = std::max(theta0, theta0_floor_rad);

    double sigma = theta0;
    if (use_detector_sigma) {
      const double segCm = xOverX0 * radiation_length_cm;
      double sigmaDet2 = 0.0;
      if (std::isfinite(segCm) && segCm > 0.0) sigmaDet2 = detector_sigma_a / (segCm * segCm) + detector_sigma_c;
      if (!std::isfinite(sigmaDet2) || sigmaDet2 <= 0.0)
        sigmaDet2 = detector_sigma_floor_rad * detector_sigma_floor_rad;
      double sigmaDet = std::sqrt(sigmaDet2);
      sigmaDet = std::max(sigmaDet, detector_sigma_floor_rad);
      sigma = std::sqrt(theta0 * theta0 + sigmaDet * sigmaDet);
    }
    sigma = std::max(sigma, theta0_floor_rad);

    const double dtxz = theta_xz[i];
    const double dtyz = theta_yz[i];
    if (!std::isfinite(dtxz) || !std::isfinite(dtyz)) return std::numeric_limits<double>::infinity();
    const double pullX = dtxz / sigma;
    const double pullY = dtyz / sigma;
    nll += 0.5 * (pullX * pullX + pullY * pullY) + 2.0 * std::log(sigma);
  }
  return (std::isfinite(nll)) ? nll : std::numeric_limits<double>::infinity();
}

double MCSLikelihood::ComputeNLL(double momentumGeV) const {
  MCSLikelihoodConfig c;
  c.radiationLengthCm = radiation_length_cm_;
  c.theta0FloorRad = theta0_floor_rad_;
  c.maxAbsDeltaThetaRad = -1.0;
  c.useDetectorSigma = use_detector_sigma_;
  c.detectorSigmaFloorRad = detector_sigma_floor_rad_;
  c.detectorSigmaA = detector_sigma_a_;
  c.detectorSigmaC = detector_sigma_c_;
  return ComputePionMcsNLLFromScatteringSamples(theta_xz_, theta_yz_, x_over_x0_, cumulative_length_cm_, momentumGeV, c);
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

bool EstimatePionMomentumMCSReco(AnaParticlePD* track, const PionMCSConfig& cfg, PionMCSResult& out,
                                 PionMCSRecoBuffers& reco, std::vector<double>* pAxisGeVOut,
                                 std::vector<double>* logLOut) {
  out.bestMomentumGeV = -999.0;
  out.bestLogLikelihood = -999.0;
  out.valid = false;
  if (pAxisGeVOut) pAxisGeVOut->clear();
  if (logLOut) logLOut->clear();
  reco.clear();
  if (!track) return false;

  std::vector<double> xoverx0;
  if (!BuildPionMcsScatteringSamples(*track, cfg.likelihood, reco.deltaThetaPair, xoverx0, nullptr, &reco.deltaThetaXz,
                                     &reco.deltaThetaYz, nullptr, nullptr, nullptr, &reco.segmentLengthSingle,
                                     &reco.segmentFitChi2NdfSingle, &reco.segmentRrToEndSingle))
    return false;

  reco.segmentLengthPair.clear();
  reco.segmentLengthPair.reserve(xoverx0.size());
  for (size_t i = 0; i < xoverx0.size(); ++i)
    reco.segmentLengthPair.push_back(xoverx0[i] * cfg.likelihood.radiationLengthCm);

  const int nMcs = static_cast<int>(reco.deltaThetaPair.size());
  int dropFirst = 0;
  int dropLast = 0;
  if (nMcs > 0) {
    dropFirst = std::max(0, cfg.dropFirstNSegments);
    dropLast = std::max(0, cfg.dropLastNSegments);
    dropFirst = std::min(dropFirst, nMcs - 1);
    dropLast = std::min(dropLast, nMcs - 1);
    if (dropFirst + dropLast >= nMcs) {
      dropFirst = 0;
      dropLast = 0;
    }
  }
  const int firstKeep = dropFirst;
  const int lastKeepExcl = nMcs - dropLast;
  if (nMcs > 0 && lastKeepExcl > firstKeep) {
    reco.scatterUsedPair.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    reco.scatterUsedXz.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    reco.scatterUsedYz.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    reco.scatterUsedRR.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    std::vector<double> xoverx0used;
    xoverx0used.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    std::vector<double> cumMidUsed;
    cumMidUsed.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
    std::vector<double> cumMidFull(static_cast<size_t>(nMcs));
    double accumMid = 0.0;
    for (int j = 0; j < nMcs; ++j) {
      const double seg = reco.segmentLengthPair[static_cast<size_t>(j)];
      cumMidFull[static_cast<size_t>(j)] = accumMid + 0.5 * seg;
      accumMid += seg;
    }
    const double totalLen = std::accumulate(reco.segmentLengthPair.begin(), reco.segmentLengthPair.end(), 0.0);
    double accumLen = 0.0;
    for (int i = 0; i < firstKeep; ++i) accumLen += reco.segmentLengthPair[static_cast<size_t>(i)];
    for (int i = firstKeep; i < lastKeepExcl; ++i) {
      const double seg = reco.segmentLengthPair[static_cast<size_t>(i)];
      reco.scatterUsedPair.push_back(reco.deltaThetaPair[static_cast<size_t>(i)]);
      reco.scatterUsedXz.push_back(reco.deltaThetaXz[static_cast<size_t>(i)]);
      reco.scatterUsedYz.push_back(reco.deltaThetaYz[static_cast<size_t>(i)]);
      const double segCenter = accumLen + 0.5 * seg;
      reco.scatterUsedRR.push_back(std::max(0.0, totalLen - segCenter));
      accumLen += seg;
      xoverx0used.push_back(xoverx0[static_cast<size_t>(i)]);
      cumMidUsed.push_back(cumMidFull[static_cast<size_t>(i)]);
    }

    std::vector<double> pAxisLocal;
    std::vector<double> logLlocal;
    std::vector<double>* pAxisPtr = pAxisGeVOut ? pAxisGeVOut : &pAxisLocal;
    std::vector<double>* logLPtr = logLOut ? logLOut : &logLlocal;
    pAxisPtr->clear();
    logLPtr->clear();
    FillMomentumScanAxis(cfg.momentumScanMaxGeV, *pAxisPtr);
    bool curveOk = true;
    for (double pGeV : *pAxisPtr) {
      if (!std::isfinite(pGeV) || pGeV <= 0.0) {
        curveOk = false;
        break;
      }
      const double nll =
          ComputePionMcsNLLFromScatteringSamples(reco.scatterUsedXz, reco.scatterUsedYz, xoverx0used, cumMidUsed, pGeV,
                                                 cfg.likelihood);
      if (!std::isfinite(nll)) {
        curveOk = false;
        break;
      }
      logLPtr->push_back(-nll);
    }
    if (curveOk && logLPtr->size() == pAxisPtr->size()) {
      double bestLogL = -std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < pAxisPtr->size(); ++i) {
        if (!std::isfinite((*logLPtr)[i])) continue;
        if ((*logLPtr)[i] > bestLogL) {
          bestLogL = (*logLPtr)[i];
          out.bestMomentumGeV = (*pAxisPtr)[i];
        }
      }
      if (std::isfinite(out.bestMomentumGeV) && out.bestMomentumGeV > 0.0) {
        out.bestLogLikelihood = std::isfinite(bestLogL) ? bestLogL : -999.0;
        out.valid = true;
      }
    }
  }

  return true;
}

bool EstimatePionMomentumMCSTrue(const std::vector<TVector3>& orderedScePositions, const PionMCSConfig& cfg,
                                 PionMCSResult& out, PionMCSTrueBuffers& buf) {
  out.bestMomentumGeV = -999.0;
  out.bestLogLikelihood = -999.0;
  out.valid = false;
  buf.clear();

  std::vector<double> maintruexoverx0;
  if (!BuildPionMcsScatteringSamplesFromOrderedPositions(orderedScePositions, cfg.likelihood, buf.deltaThetaPair,
                                                        maintruexoverx0, nullptr, &buf.deltaThetaXz, &buf.deltaThetaYz,
                                                        nullptr, nullptr, nullptr, &buf.segmentLengthSingle,
                                                        &buf.segmentFitChi2NdfSingle, &buf.segmentRrToEndSingle))
    return false;

  buf.segmentLengthPair.clear();
  buf.segmentLengthPair.reserve(maintruexoverx0.size());
  for (size_t i = 0; i < maintruexoverx0.size(); ++i)
    buf.segmentLengthPair.push_back(maintruexoverx0[i] * cfg.likelihood.radiationLengthCm);

  const int nMcsT = static_cast<int>(buf.deltaThetaPair.size());
  int dropFirstT = 0;
  int dropLastT = 0;
  if (nMcsT > 0) {
    dropFirstT = std::max(0, cfg.dropFirstNSegments);
    dropLastT = std::max(0, cfg.dropLastNSegments);
    dropFirstT = std::min(dropFirstT, nMcsT - 1);
    dropLastT = std::min(dropLastT, nMcsT - 1);
    if (dropFirstT + dropLastT >= nMcsT) {
      dropFirstT = 0;
      dropLastT = 0;
    }
  }
  const int firstKeepT = dropFirstT;
  const int lastKeepExclT = nMcsT - dropLastT;
  if (nMcsT > 0 && lastKeepExclT > firstKeepT) {
    std::vector<double> xzusd;
    std::vector<double> yzusd;
    std::vector<double> xoverx0usedT;
    std::vector<double> cumMidUsedT;
    xzusd.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
    yzusd.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
    xoverx0usedT.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
    cumMidUsedT.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
    std::vector<double> cumMidFullT(static_cast<size_t>(nMcsT));
    double accMidT = 0.0;
    for (int j = 0; j < nMcsT; ++j) {
      const double seg = buf.segmentLengthPair[static_cast<size_t>(j)];
      cumMidFullT[static_cast<size_t>(j)] = accMidT + 0.5 * seg;
      accMidT += seg;
    }
    for (int i = firstKeepT; i < lastKeepExclT; ++i) {
      xzusd.push_back(buf.deltaThetaXz[static_cast<size_t>(i)]);
      yzusd.push_back(buf.deltaThetaYz[static_cast<size_t>(i)]);
      xoverx0usedT.push_back(maintruexoverx0[static_cast<size_t>(i)]);
      cumMidUsedT.push_back(cumMidFullT[static_cast<size_t>(i)]);
    }
    std::vector<double> pAxisMcsTrueGeV;
    std::vector<double> logLMcsTrue;
    FillMomentumScanAxis(cfg.momentumScanMaxGeV, pAxisMcsTrueGeV);
    bool curveOkT = true;
    for (double pGeV : pAxisMcsTrueGeV) {
      if (!std::isfinite(pGeV) || pGeV <= 0.0) {
        curveOkT = false;
        break;
      }
      const double nll =
          ComputePionMcsNLLFromScatteringSamples(xzusd, yzusd, xoverx0usedT, cumMidUsedT, pGeV, cfg.likelihood);
      if (!std::isfinite(nll)) {
        curveOkT = false;
        break;
      }
      logLMcsTrue.push_back(-nll);
    }
    if (curveOkT && logLMcsTrue.size() == pAxisMcsTrueGeV.size()) {
      double bestLogLT = -std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < pAxisMcsTrueGeV.size(); ++i) {
        if (!std::isfinite(logLMcsTrue[i])) continue;
        if (logLMcsTrue[i] > bestLogLT) {
          bestLogLT = logLMcsTrue[i];
          out.bestMomentumGeV = pAxisMcsTrueGeV[i];
        }
      }
      if (std::isfinite(out.bestMomentumGeV) && out.bestMomentumGeV > 0.0) {
        out.bestLogLikelihood = std::isfinite(bestLogLT) ? bestLogLT : -999.0;
        out.valid = true;
      }
    }
  }

  return true;
}

}  // namespace pdMomReconstruction
