#include "pdAnnihilationUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdJointK0sPionMomentum.hxx"
#include "pdMomReconstruction.hxx"
#include "pdUtilsDEdx.hxx"
#include "neutralKaonAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

namespace pdAnnihilationUtils {
namespace {

bool HasValidStartPosition(const AnaParticlePD* particle) {
  if (!particle) return false;
  return (particle->PositionStart[0] > -900.0 &&
          particle->PositionStart[1] > -900.0 &&
          particle->PositionStart[2] > -900.0);
}

bool HasValidEndPosition(const AnaParticlePD* particle) {
  if (!particle) return false;
  return (particle->PositionEnd[0] > -900.0 &&
          particle->PositionEnd[1] > -900.0 &&
          particle->PositionEnd[2] > -900.0);
}

bool HasValidPosition3(const Float_t pos[3]) {
  if (!pos) return false;
  return (pos[0] > -900.0f && pos[1] > -900.0f && pos[2] > -900.0f);
}

// Endpoint pairing: combo 0=SS, 1=SE, 2=ES, 3=EE; combo -1 if no valid pair.
struct AnnihilationEndpointPairing {
  float distance = std::numeric_limits<float>::infinity();
  int combo = -1;
};

AnnihilationEndpointPairing BestAnnihilationPairEndpoint(const AnaParticlePD* d1, const AnaParticlePD* d2) {
  AnnihilationEndpointPairing out;
  if (!d1 || !d2) return out;
  const bool v1s = HasValidStartPosition(d1);
  const bool v1e = HasValidEndPosition(d1);
  const bool v2s = HasValidStartPosition(d2);
  const bool v2e = HasValidEndPosition(d2);

  auto consider = [&](int comboId, const Float_t a[4], const Float_t b[4]) {
    if (!HasValidPosition3(a) || !HasValidPosition3(b)) return;
    const TVector3 va(a[0], a[1], a[2]);
    const TVector3 vb(b[0], b[1], b[2]);
    const float dist = static_cast<float>((va - vb).Mag());
    if (dist < out.distance) {
      out.distance = dist;
      out.combo = comboId;
    }
  };

  if (v1s && v2s) consider(0, d1->PositionStart, d2->PositionStart);
  if (v1s && v2e) consider(1, d1->PositionStart, d2->PositionEnd);
  if (v1e && v2s) consider(2, d1->PositionEnd, d2->PositionStart);
  if (v1e && v2e) consider(3, d1->PositionEnd, d2->PositionEnd);

  return out;
}

// Start–start distance only (legacy/plain pairing). combo=SS so no ApplyAnnihilationPairingReversal.
AnnihilationEndpointPairing StartStartAnnihilationPairing(const AnaParticlePD* d1, const AnaParticlePD* d2) {
  AnnihilationEndpointPairing out;
  if (!d1 || !d2) return out;
  if (!HasValidStartPosition(d1) || !HasValidStartPosition(d2)) return out;
  const TVector3 a(d1->PositionStart[0], d1->PositionStart[1], d1->PositionStart[2]);
  const TVector3 b(d2->PositionStart[0], d2->PositionStart[1], d2->PositionStart[2]);
  out.distance = static_cast<float>((a - b).Mag());
  out.combo = 0;
  return out;
}

// Swap start/end, reverse trajectory and hit order, flip TrajectoryDirection when start Z > end Z.
void ReverseParticleIfStartZGreaterThanEndZ(AnaParticlePD* p) {
  if (!p) return;
  if (!HasValidStartPosition(p) || !HasValidEndPosition(p)) return;
  if (p->PositionStart[2] <= p->PositionEnd[2]) return;

  for (int k = 0; k < 4; ++k) {
    std::swap(p->PositionStart[k], p->PositionEnd[k]);
  }
  for (int k = 0; k < 3; ++k) {
    std::swap(p->DirectionStart[k], p->DirectionEnd[k]);
  }
  for (int k = 0; k < 4; ++k) {
    std::swap(p->PositionStartSCE[k], p->PositionEndSCE[k]);
  }
  for (int k = 0; k < 3; ++k) {
    std::swap(p->DirectionStartSCE[k], p->DirectionEndSCE[k]);
  }
  std::swap(p->Momentum, p->MomentumEnd);

  std::reverse(p->TrjPoints.begin(), p->TrjPoints.end());
  for (int pl = 0; pl < 3; ++pl) {
    std::reverse(p->Hits[pl].begin(), p->Hits[pl].end());
  }
  if (p->TrajectoryDirectionNPoints > 0) {
    p->TrajectoryDirection = -1.0 * p->TrajectoryDirection;
  }
  if (p->Hits[2].size() >= 2) {
    pdAnaUtils::EstimateHitsDirection(p);
    pdAnaUtils::ComputeResidualRange(p);
  }
}

// SS: no reversal. SE/ES/EE: reverse each daughter if PositionStart(Z) > PositionEnd(Z).
void ApplyAnnihilationPairingReversal(int combo, AnaParticlePD* d1, AnaParticlePD* d2) {
  if (combo <= 0) return;
  ReverseParticleIfStartZGreaterThanEndZ(d1);
  ReverseParticleIfStartZGreaterThanEndZ(d2);
}

TVector3 GetAnnihilationVertexPositionForDegeneracy(const AnaAnnihilationVertexPD* vertex) {
  if (vertex && HasValidPosition3(vertex->PositionFit)) {
    return TVector3(vertex->PositionFit[0], vertex->PositionFit[1], vertex->PositionFit[2]);
  }
  return TVector3(-999.0, -999.0, -999.0);
}

TVector3 GetCreationVertexPositionForDegeneracy(const AnaCreationVertexPD* vertex) {
  if (vertex && HasValidPosition3(vertex->Position)) {
    return TVector3(vertex->Position[0], vertex->Position[1], vertex->Position[2]);
  }
  if (vertex && HasValidPosition3(vertex->PositionPandora)) {
    return TVector3(vertex->PositionPandora[0], vertex->PositionPandora[1], vertex->PositionPandora[2]);
  }
  return TVector3(-999.0, -999.0, -999.0);
}

bool HasValidPoint(const TVector3& point) {
  return std::isfinite(point.X()) && std::isfinite(point.Y()) && std::isfinite(point.Z()) &&
         point.X() > -900.0 && point.Y() > -900.0 && point.Z() > -900.0;
}

TVector3 ProjectPointOntoLine(const TVector3& point, const TVector3& linePoint, const TVector3& lineDirection) {
  if (!HasValidPoint(point) || !HasValidPoint(linePoint)) return TVector3(-999.0, -999.0, -999.0);
  TVector3 dir = lineDirection;
  if (!std::isfinite(dir.X()) || !std::isfinite(dir.Y()) || !std::isfinite(dir.Z()) || dir.Mag2() <= 1e-10) {
    return TVector3(-999.0, -999.0, -999.0);
  }
  dir = dir.Unit();
  return linePoint + (point - linePoint).Dot(dir) * dir;
}

bool GetParticleFitLine(const AnaParticlePD* particle,
                        double trackFitLength,
                        double trackFitDistanceFromStart,
                        TVector3& fitAnchor,
                        TVector3& fitDir) {
  std::vector<double> fitParams;
  pdAnaUtils::ExtrapolateTrack(const_cast<AnaParticlePD*>(particle), fitParams, trackFitLength, true,
                               trackFitDistanceFromStart);
  const bool fitValid = (fitParams.size() >= 6 && std::isfinite(fitParams[0]) && std::isfinite(fitParams[1]) &&
                         std::isfinite(fitParams[2]) && std::isfinite(fitParams[3]) && std::isfinite(fitParams[4]) &&
                         std::isfinite(fitParams[5]) && fitParams[0] > -900.0 && fitParams[1] > -900.0 &&
                         fitParams[2] > -900.0);
  if (!fitValid) return false;
  fitAnchor.SetXYZ(fitParams[0], fitParams[1], fitParams[2]);
  fitDir.SetXYZ(fitParams[3], fitParams[4], fitParams[5]);
  return HasValidPoint(fitAnchor) && fitDir.Mag2() > 1e-10;
}

double EstimatePathDistanceFromStart(const AnaParticlePD* particle, const TVector3& position) {
  if (!particle || !HasValidPoint(position)) return -1.0;

  std::vector<std::pair<TVector3, double>> trajectoryPointsWithDistance;
  if (particle->TrjPoints.size() >= 2) {
    trajectoryPointsWithDistance.reserve(particle->TrjPoints.size());
    double cumulative = 0.0;
    TVector3 prev;
    bool hasPrev = false;
    for (size_t i = 0; i < particle->TrjPoints.size(); ++i) {
      const TVector3 pos = particle->TrjPoints[i].Position;
      if (!HasValidPoint(pos)) continue;
      if (hasPrev) cumulative += (pos - prev).Mag();
      trajectoryPointsWithDistance.push_back(std::make_pair(pos, cumulative));
      prev = pos;
      hasPrev = true;
    }
  }

  if (!trajectoryPointsWithDistance.empty()) {
    double bestDist2 = 1e30;
    double bestArc = -1.0;
    for (const auto& tp : trajectoryPointsWithDistance) {
      const double d2 = (position - tp.first).Mag2();
      if (d2 < bestDist2) {
        bestDist2 = d2;
        bestArc = tp.second;
      }
    }
    return bestArc;
  }

  if (!HasValidStartPosition(particle)) return -1.0;
  const TVector3 referencePos(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
  TVector3 travelDir(particle->DirectionStart[0], particle->DirectionStart[1], particle->DirectionStart[2]);
  if (travelDir.Mag2() > 1e-10) {
    travelDir = travelDir.Unit();
    return (position - referencePos).Dot(travelDir);
  }
  return (position - referencePos).Mag();
}

bool ParticleHasRawTailSupportNearVertex(const AnaParticlePD* particle,
                                         const TVector3& vertexPos,
                                         double vertexRadius,
                                         double trackFitDistanceFromStart) {
  if (!particle || !HasValidPoint(vertexPos) || vertexRadius <= 0.0) return false;
  if (HasValidStartPosition(particle)) {
    const TVector3 rawStart(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
    if ((rawStart - vertexPos).Mag() <= vertexRadius) return true;
  }
  for (const AnaHitPD& hit : particle->Hits[2]) {
    const TVector3 hitPos = hit.Position;
    if (!HasValidPoint(hitPos)) continue;
    const double pathDistance = EstimatePathDistanceFromStart(particle, hitPos);
    if (pathDistance < 0.0 || pathDistance > trackFitDistanceFromStart) continue;
    if ((hitPos - vertexPos).Mag() <= vertexRadius) return true;
  }
  return false;
}

bool ParticleHasProjectedTailSupportNearPoint(const AnaParticlePD* particle,
                                              const TVector3& fitAnchor,
                                              const TVector3& fitDir,
                                              const TVector3& referencePoint,
                                              double maxDistance,
                                              double trackFitDistanceFromStart) {
  if (!particle || !HasValidPoint(fitAnchor) || fitDir.Mag2() <= 1e-10 || !HasValidPoint(referencePoint) ||
      maxDistance <= 0.0) {
    return false;
  }
  if (HasValidStartPosition(particle)) {
    const TVector3 rawStart(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
    const TVector3 projectedStart = ProjectPointOntoLine(rawStart, fitAnchor, fitDir);
    if (HasValidPoint(projectedStart) && (projectedStart - referencePoint).Mag() <= maxDistance) return true;
  }
  for (const AnaHitPD& hit : particle->Hits[2]) {
    const TVector3 hitPos = hit.Position;
    if (!HasValidPoint(hitPos)) continue;
    const double pathDistance = EstimatePathDistanceFromStart(particle, hitPos);
    if (pathDistance < 0.0 || pathDistance > trackFitDistanceFromStart) continue;
    const TVector3 projectedHit = ProjectPointOntoLine(hitPos, fitAnchor, fitDir);
    if (HasValidPoint(projectedHit) && (projectedHit - referencePoint).Mag() <= maxDistance) return true;
  }
  return false;
}

bool IsVertexDaughter(const AnaAnnihilationVertexPD* vertex, const AnaParticlePD* candidate) {
  if (!vertex || !candidate) return false;
  for (AnaParticlePD* daughter : vertex->Particles) {
    if (daughter == candidate) return true;
  }
  return false;
}

Int_t ComputeVertexDegeneracyAtPosition(const AnaEventB& event,
                                        const TVector3& referencePos,
                                        double vertexRadius,
                                        double lineToVertexDistance,
                                        double originSupportDistance,
                                        double trackFitLength,
                                        double trackFitDistanceFromStart,
                                        const std::unordered_set<Int_t>& excludedUniqueIds,
                                        const AnaAnnihilationVertexPD* daughterExclusionVertex = nullptr,
                                        const TVector3* annihilationPosForVeto = nullptr,
                                        const TVector3* creationPosForVeto = nullptr) {
  if (!HasValidPoint(referencePos) || vertexRadius <= 0.0 || lineToVertexDistance <= 0.0 ||
      originSupportDistance <= 0.0) {
    return 0;
  }

  const bool applyCreationVsAnnihilationVeto =
      annihilationPosForVeto && creationPosForVeto && HasValidPoint(*annihilationPosForVeto) &&
      HasValidPoint(*creationPosForVeto);

  Int_t degeneracy = 0;
  for (Int_t p = 0; p < event.nParticles; ++p) {
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[p]);
    if (!particle) continue;
    if (!excludedUniqueIds.empty() && excludedUniqueIds.count(particle->UniqueID) > 0) continue;
    if (daughterExclusionVertex && IsVertexDaughter(daughterExclusionVertex, particle)) continue;
    if (!ParticleHasRawTailSupportNearVertex(particle, referencePos, vertexRadius, trackFitDistanceFromStart)) {
      continue;
    }

    TVector3 fitAnchor;
    TVector3 fitDir;
    if (!GetParticleFitLine(particle, trackFitLength, trackFitDistanceFromStart, fitAnchor, fitDir)) continue;
    fitDir = fitDir.Unit();

    const TVector3 closestPointToReference = ProjectPointOntoLine(referencePos, fitAnchor, fitDir);
    if (!HasValidPoint(closestPointToReference)) continue;
    if ((closestPointToReference - referencePos).Mag() > lineToVertexDistance) continue;

    if (!ParticleHasProjectedTailSupportNearPoint(particle, fitAnchor, fitDir, closestPointToReference,
                                                  originSupportDistance, trackFitDistanceFromStart)) {
      continue;
    }

    if (applyCreationVsAnnihilationVeto) {
      const TVector3 closestToAnnihilation = ProjectPointOntoLine(*annihilationPosForVeto, fitAnchor, fitDir);
      const TVector3 closestToCreation = ProjectPointOntoLine(*creationPosForVeto, fitAnchor, fitDir);
      if (HasValidPoint(closestToAnnihilation) && HasValidPoint(closestToCreation)) {
        const double dAnnihilation = (closestToAnnihilation - *annihilationPosForVeto).Mag();
        const double dCreation = (closestToCreation - *creationPosForVeto).Mag();
        if (dCreation < dAnnihilation) continue;
      }
    }

    ++degeneracy;
  }

  return degeneracy;
}

bool IsProtonLikeAndNotPionLike(const AnaParticlePD* particle, double protonChi2NdfRejectBelow,
                                double pionChi2NdfRejectAbove) {
  if (!particle) return false;

  const Float_t protonChi2Ndf = pdAnaUtils::Chi2PIDChi2PerHit(particle, 2212);
  if (!(protonChi2Ndf > -900.f)) return false;

  if (static_cast<double>(protonChi2Ndf) < protonChi2NdfRejectBelow) return true;

  const Float_t pionChi2Ndf = pdAnaUtils::Chi2PIDChi2PerHit(particle, 211);
  if (!(pionChi2Ndf > -900.f)) return false;

  return static_cast<double>(pionChi2Ndf) > pionChi2NdfRejectAbove;
}

enum DaughterMomentumMethod {
  kMomMethodUnassigned = -1,
  kMomMethodPionRangeStopping = 0,
  kMomMethodCalorimetric = 1,
  kMomMethodFailed = 2,
  kMomMethodFreeRangeML = 3,
  kMomMethodJointK0sGrid = 4
};

struct DaughterMomentumDebugInfo {
  Int_t hasPreexistingMomentum = -1;
  Int_t extensionAttempted = 0;
  Int_t extensionValid = 0;
  Float_t extensionChi2Ndf = -999.0f;
  Int_t extensionNValidHits = 0;
  Float_t extensionDedxBias = -999.0f;   // Gaussian mean of (measured - Bethe-Bloch) dEdx [MeV/cm]
  Float_t extensionDedxSigma = -999.0f;  // Gaussian sigma of (measured - Bethe-Bloch) dEdx [MeV/cm]
  Int_t extensionDedxFitOk = -1;
};

Int_t CountValidCollectionPlaneHits(const AnaParticlePD* particle) {
  if (!particle) return 0;
  Int_t count = 0;
  for (size_t i = 0; i < particle->Hits[2].size(); ++i) {
    const AnaHitPD& hit = particle->Hits[2][i];
    if (hit.dEdx > 0.0f && hit.dEdx != -999.0f && hit.dEdx < 1000.0f) {
      ++count;
    }
  }
  return count;
}

Int_t AssignMomentumFromResidualRange(AnaParticlePD* particle, Int_t pdgHypothesis,
                                      DaughterMomentumDebugInfo* debugInfo = nullptr) {
  if (!particle) return kMomMethodUnassigned;

  const double scanLmaxCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                                : 450.;
  const double scanStepCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepCm")
                                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepCm")
                                : 1.;
  const int minInteriorHits =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
          ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
          : 15;
  const int skipFirst =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
          ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
          : 3;
  const int skipLast =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
          ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
          : 3;
  double dedxMinMeVcm =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
          ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
          : 0.5;
  double dedxMaxMeVcm =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
          ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
          : 5.0;
  if (pdgHypothesis != 211) {
    dedxMinMeVcm = 0.;
    dedxMaxMeVcm = 0.;
  }
  const double pdfPathCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                               ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                               : 0.65;
  const double scanStepFineCm =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepFineCm")
          ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepFineCm")
          : 0.;
  const double lowPMomentumRefineGeV =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
          ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
          : 0.2;
  const bool computeDedxBiasDiagnostics =
      !ND::params().HasParameter("neutralKaonAnalysis.FreeRangeComputeDedxBiasDiagnostics") ||
      ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeComputeDedxBiasDiagnostics") != 0;

  const Int_t nCollHits = static_cast<Int_t>(particle->Hits[2].size());
  const Int_t minHitsOnTrack = skipFirst + skipLast + std::max(1, minInteriorHits);
  const Int_t attempted = (nCollHits >= minHitsOnTrack) ? 1 : 0;
  if (debugInfo) {
    debugInfo->hasPreexistingMomentum = -1;
    debugInfo->extensionAttempted = attempted;
    debugInfo->extensionValid = 0;
    debugInfo->extensionChi2Ndf = -999.0f;
    debugInfo->extensionDedxBias = -999.0f;
    debugInfo->extensionDedxSigma = -999.0f;
    debugInfo->extensionDedxFitOk = -1;
    debugInfo->extensionNValidHits = CountValidCollectionPlaneHits(particle);
  }

  const pdAnaUtils::DEdxFreeRangeFitResult fit = pdAnaUtils::GetdEdxLikelihoodFreeRangeFit(
      particle, pdgHypothesis, scanLmaxCm, scanStepCm, minInteriorHits, skipFirst, skipLast, dedxMinMeVcm, dedxMaxMeVcm,
      pdfPathCm, computeDedxBiasDiagnostics, scanStepFineCm, lowPMomentumRefineGeV);

  if (debugInfo && std::isfinite(static_cast<double>(fit.logLikelihood))) {
    debugInfo->extensionChi2Ndf = static_cast<Float_t>(fit.logLikelihood);
  }
  if (debugInfo && std::isfinite(static_cast<double>(fit.meanDedxBias))) {
    debugInfo->extensionDedxBias = static_cast<Float_t>(fit.meanDedxBias);
  }
  if (debugInfo && std::isfinite(static_cast<double>(fit.sigmaDedxBias))) {
    debugInfo->extensionDedxSigma = static_cast<Float_t>(fit.sigmaDedxBias);
  }
  if (debugInfo) {
    debugInfo->extensionDedxFitOk = fit.dedxFitOk;
  }

  if (std::isfinite(static_cast<double>(fit.momentumGeV)) && fit.momentumGeV > 0.f) {
    particle->Momentum = fit.momentumGeV;
    if (debugInfo) debugInfo->extensionValid = 1;
    return kMomMethodFreeRangeML;
  }

  if (debugInfo) debugInfo->extensionValid = 0;
  return kMomMethodFailed;
}

Int_t AssignPionMomentumFromResidualRange(AnaParticlePD* particle, DaughterMomentumDebugInfo* debugInfo = nullptr) {
  return AssignMomentumFromResidualRange(particle, 211, debugInfo);
}

void DaughterPandoraAndFitDirs(AnaAnnihilationVertexPD* vertex, double trackFitLength, double trackFitDistanceFromStart,
                               TVector3& dirPandora1, TVector3& dirPandora2, TVector3& dirFit1, TVector3& dirFit2) {
  AnaParticlePD* daughter1 = vertex->Particles[0];
  AnaParticlePD* daughter2 = vertex->Particles[1];
  dirPandora1.SetXYZ(daughter1->DirectionStart[0], daughter1->DirectionStart[1], daughter1->DirectionStart[2]);
  dirPandora2.SetXYZ(daughter2->DirectionStart[0], daughter2->DirectionStart[1], daughter2->DirectionStart[2]);
  dirFit1 = dirPandora1;
  dirFit2 = dirPandora2;
  std::vector<double> fit1;
  std::vector<double> fit2;
  pdAnaUtils::ExtrapolateTrack(daughter1, fit1, trackFitLength, true, trackFitDistanceFromStart);
  pdAnaUtils::ExtrapolateTrack(daughter2, fit2, trackFitLength, true, trackFitDistanceFromStart);
  const bool fit1Valid = (fit1.size() >= 6 && fit1[3] > -900.0 && fit1[4] > -900.0 && fit1[5] > -900.0 &&
                          std::isfinite(fit1[3]) && std::isfinite(fit1[4]) && std::isfinite(fit1[5]));
  const bool fit2Valid = (fit2.size() >= 6 && fit2[3] > -900.0 && fit2[4] > -900.0 && fit2[5] > -900.0 &&
                          std::isfinite(fit2[3]) && std::isfinite(fit2[4]) && std::isfinite(fit2[5]));
  if (fit1Valid) dirFit1.SetXYZ(fit1[3], fit1[4], fit1[5]);
  if (fit2Valid) dirFit2.SetXYZ(fit2[3], fit2[4], fit2[5]);
}

bool ComputeVertexKinematicsWithDirections(Float_t p1Mag, Float_t p2Mag, TVector3 dir1, TVector3 dir2,
                                           Float_t& momentumOut, Float_t& invariantMassOut) {
  static constexpr Float_t kPionMassGeV = 0.13957f;
  if (!(std::isfinite(p1Mag) && std::isfinite(p2Mag) && p1Mag > 0.0f && p2Mag > 0.0f)) return false;
  if (dir1.Mag2() <= 0.0 || dir2.Mag2() <= 0.0) return false;

  dir1 = dir1.Unit();
  dir2 = dir2.Unit();
  const TVector3 p1 = p1Mag * dir1;
  const TVector3 p2 = p2Mag * dir2;
  const TVector3 pTot = p1 + p2;
  const double e1 = std::sqrt(static_cast<double>(p1Mag) * p1Mag +
                              static_cast<double>(kPionMassGeV) * kPionMassGeV);
  const double e2 = std::sqrt(static_cast<double>(p2Mag) * p2Mag +
                              static_cast<double>(kPionMassGeV) * kPionMassGeV);
  const double eTot = e1 + e2;
  const double m2 = eTot * eTot - pTot.Mag2();

  momentumOut = static_cast<Float_t>(pTot.Mag());
  invariantMassOut = (m2 > 0.0) ? static_cast<Float_t>(std::sqrt(m2)) : 0.0f;
  return true;
}

void FillVertexKinematicsFromDaughters(AnaAnnihilationVertexPD* vertex, double trackFitLength,
                                       double trackFitDistanceFromStart) {
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* daughter1 = vertex->Particles[0];
  AnaParticlePD* daughter2 = vertex->Particles[1];
  if (!daughter1 || !daughter2) return;

  const Float_t p1Mag = daughter1->Momentum;
  const Float_t p2Mag = daughter2->Momentum;

  TVector3 dirPandora1, dirPandora2, dirFit1, dirFit2;
  DaughterPandoraAndFitDirs(vertex, trackFitLength, trackFitDistanceFromStart, dirPandora1, dirPandora2, dirFit1,
                            dirFit2);
  ComputeVertexKinematicsWithDirections(p1Mag, p2Mag, dirPandora1, dirPandora2,
                                        vertex->MomentumPandora, vertex->InvariantMassPandora);
  ComputeVertexKinematicsWithDirections(p1Mag, p2Mag, dirFit1, dirFit2,
                                        vertex->MomentumFit, vertex->InvariantMassFit);
  // Preserve previous API behavior using the fit-based values.
  vertex->Momentum = vertex->MomentumFit;
  vertex->InvariantMass = vertex->InvariantMassFit;
}

Float_t AngleBetweenUnitVectors(const TVector3& a, const TVector3& b) {
  const double c = std::max(-1.0, std::min(1.0, a.Unit().Dot(b.Unit())));
  return static_cast<Float_t>(std::acos(c));
}

Int_t ComputeAnnihilationVertexDegeneracy(const AnaEventB& event,
                                          const AnaAnnihilationVertexPD* vertex,
                                          double annihilationVertexRadius,
                                          double annihilationVertexLineToVertexDistance,
                                          double annihilationVertexOriginSupportDistance,
                                          double trackFitLength,
                                          double trackFitDistanceFromStart,
                                          Int_t excludedParticleUniqueID = -1,
                                          const AnaCreationVertexPD* creationVertex = nullptr) {
  if (!vertex || annihilationVertexRadius <= 0.0 || annihilationVertexLineToVertexDistance <= 0.0 ||
      annihilationVertexOriginSupportDistance <= 0.0) {
    return 0;
  }

  const TVector3 annihilationPos = GetAnnihilationVertexPositionForDegeneracy(vertex);
  if (!HasValidPoint(annihilationPos)) return 0;

  std::unordered_set<Int_t> excludedUniqueIds;
  if (excludedParticleUniqueID >= 0) {
    excludedUniqueIds.insert(excludedParticleUniqueID);
  }
  const TVector3 creationPos = GetCreationVertexPositionForDegeneracy(creationVertex);

  return ComputeVertexDegeneracyAtPosition(event, annihilationPos, annihilationVertexRadius,
                                           annihilationVertexLineToVertexDistance,
                                           annihilationVertexOriginSupportDistance, trackFitLength,
                                           trackFitDistanceFromStart, excludedUniqueIds, vertex,
                                           &annihilationPos, HasValidPoint(creationPos) ? &creationPos : nullptr);
}

Int_t ComputeCreationVertexDegeneracy(const AnaEventB& event,
                                      const AnaCreationVertexPD* creationVertex,
                                      const AnaAnnihilationVertexPD* annihilationVertex,
                                      double creationVertexRadius,
                                      double creationVertexLineToVertexDistance,
                                      double creationVertexOriginSupportDistance,
                                      double trackFitLength,
                                      double trackFitDistanceFromStart,
                                      Int_t excludedParticleUniqueID = -1) {
  if (!creationVertex || creationVertexRadius <= 0.0 || creationVertexLineToVertexDistance <= 0.0 ||
      creationVertexOriginSupportDistance <= 0.0) {
    return 0;
  }
  const TVector3 creationPos = GetCreationVertexPositionForDegeneracy(creationVertex);
  if (!HasValidPoint(creationPos)) return 0;

  std::unordered_set<Int_t> excludedUniqueIds;
  if (excludedParticleUniqueID >= 0) {
    excludedUniqueIds.insert(excludedParticleUniqueID);
  }
  if (annihilationVertex) {
    for (AnaParticlePD* daughter : annihilationVertex->Particles) {
      if (daughter) excludedUniqueIds.insert(daughter->UniqueID);
    }
  }

  return ComputeVertexDegeneracyAtPosition(event, creationPos, creationVertexRadius,
                                           creationVertexLineToVertexDistance,
                                           creationVertexOriginSupportDistance, trackFitLength,
                                           trackFitDistanceFromStart, excludedUniqueIds, nullptr,
                                           nullptr, nullptr);
}

} // namespace

namespace {

constexpr double kJointK0sPionMassMeV = 139.57;

bool JointK0sCollectionResidualRangeLooksUnset(const AnaParticlePD* part) {
  if (!part || part->Hits[2].empty()) return false;
  for (const AnaHitPD& h : part->Hits[2]) {
    if (h.ResidualRange > 0 && std::isfinite(static_cast<double>(h.ResidualRange))) return false;
  }
  return true;
}

bool JointK0sInteriorDedxRr(AnaParticlePD* part, int minInteriorHits, int skipFirst, int skipLast,
                           double dedxMinMeVcm, double dedxMaxMeVcm, std::vector<double>& dedx,
                           std::vector<double>& rr) {
  dedx.clear();
  rr.clear();
  if (!part || part->Hits[2].empty()) return false;
  if (skipFirst < 0) skipFirst = 0;
  if (skipLast < 0) skipLast = 0;
  const int n = static_cast<int>(part->Hits[2].size());
  if (n < skipFirst + skipLast + 1) return false;
  for (int ihit = skipFirst; ihit < n - skipLast; ++ihit) {
    const AnaHitPD& h = part->Hits[2][ihit];
    dedx.push_back(static_cast<double>(h.dEdx));
    rr.push_back(static_cast<double>(h.ResidualRange));
  }
  const bool dedxWindow = (dedxMinMeVcm > 0. && dedxMaxMeVcm > dedxMinMeVcm);
  if (dedxWindow) {
    std::vector<double> nd, nr;
    nd.reserve(dedx.size());
    nr.reserve(rr.size());
    for (size_t i = 0; i < dedx.size(); ++i) {
      if (dedx[i] >= dedxMinMeVcm && dedx[i] <= dedxMaxMeVcm) {
        nd.push_back(dedx[i]);
        nr.push_back(rr[i]);
      }
    }
    dedx.swap(nd);
    rr.swap(nr);
  }
  return static_cast<int>(dedx.size()) >= minInteriorHits;
}

double JointK0sMeasuredLenCm(const AnaParticlePD* part, const std::vector<double>& rrInterior) {
  if (part && part->Length > 0.f && part->Length != -999.f) return static_cast<double>(part->Length);
  double mx = 0.;
  for (double r : rrInterior) {
    if (r > mx) mx = r;
  }
  return mx;
}

bool JointK0sPrepPionInterior(AnaParticlePD* part, int minInteriorHits, int skipFirst, int skipLast,
                             double dedxMinMeVcm, double dedxMaxMeVcm, std::vector<double>& dedx,
                             std::vector<double>& rr) {
  if (JointK0sCollectionResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);
  return JointK0sInteriorDedxRr(part, minInteriorHits, skipFirst, skipLast, dedxMinMeVcm, dedxMaxMeVcm, dedx, rr);
}

bool JointK0sBestExtensionForMomentum(const AnaParticlePD* part, const std::vector<double>& rrInterior,
                                     double pTargetGeV, double scanLmaxCm, TGraph* ke, double& extOut) {
  extOut = 0.;
  if (!part || !ke || !(pTargetGeV > 0.) || !std::isfinite(pTargetGeV) || !(scanLmaxCm > 0.)) return false;
  const double len = JointK0sMeasuredLenCm(part, rrInterior);
  double bestDiff = 1e300;
  double bestExt = 0.;
  constexpr int kNScan = 128;
  for (int k = 0; k <= kNScan; ++k) {
    const double ext = scanLmaxCm * (static_cast<double>(k) / static_cast<double>(kNScan));
    const double p = pdMomShared::RangeCmToMomentumGeV(len + ext, 211, ke, kJointK0sPionMassMeV);
    if (!std::isfinite(p)) continue;
    const double diff = std::abs(p - pTargetGeV);
    if (diff < bestDiff) {
      bestDiff = diff;
      bestExt = ext;
    }
  }
  if (!(bestDiff < 1e299)) return false;
  extOut = bestExt;
  return true;
}

} // namespace

//***************************************************************
Int_t AssignProtonMomentumFromResidualRange(AnaParticlePD* particle) {
//***************************************************************
  return AssignMomentumFromResidualRange(particle, 2212, nullptr);
}

//***************************************************************
Int_t ComputeAnnihilationVertexDegeneracyWithExclusion(const AnaEventB& event,
                                                       const AnaAnnihilationVertexPD* vertex,
                                                       Int_t excludedParticleUniqueID,
                                                       const AnaCreationVertexPD* creationVertex) {
//***************************************************************
  const double annihilationDegeneracyRadius =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          : 0.0;
  const double annihilationDegeneracyLineToVertexDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          : 0.0;
  const double annihilationDegeneracyOriginSupportDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          : 0.5 * annihilationDegeneracyLineToVertexDistance;
  const double trackFitLength =
      ND::params().HasParameter("neutralKaonAnalysis.TrackFitLength")
          ? ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength")
          : 0.0;
  const double trackFitDistanceFromStart =
      ND::params().HasParameter("neutralKaonAnalysis.TrackFitDistanceFromStart")
          ? ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart")
          : 0.0;

  return ComputeAnnihilationVertexDegeneracy(
      event, vertex, annihilationDegeneracyRadius, annihilationDegeneracyLineToVertexDistance,
      annihilationDegeneracyOriginSupportDistance, trackFitLength, trackFitDistanceFromStart,
      excludedParticleUniqueID, creationVertex);
}

//***************************************************************
Int_t ComputeCreationVertexDegeneracy(const AnaEventB& event,
                                      const AnaCreationVertexPD* creationVertex,
                                      const AnaAnnihilationVertexPD* annihilationVertex,
                                      Int_t excludedParticleUniqueID) {
//***************************************************************
  const double degeneracyRadius =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          : 0.0;
  const double lineToVertexDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          : 0.0;
  const double originSupportDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          : 0.5 * lineToVertexDistance;
  const double trackFitLength =
      ND::params().HasParameter("neutralKaonAnalysis.TrackFitLength")
          ? ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength")
          : 0.0;
  const double trackFitDistanceFromStart =
      ND::params().HasParameter("neutralKaonAnalysis.TrackFitDistanceFromStart")
          ? ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart")
          : 0.0;

  return ComputeCreationVertexDegeneracy(event, creationVertex, annihilationVertex, degeneracyRadius,
                                         lineToVertexDistance, originSupportDistance, trackFitLength,
                                         trackFitDistanceFromStart, excludedParticleUniqueID);
}

//***************************************************************
void FillNeutralParticleAlignment(AnaNeutralParticlePD* neutral, const AnaEventB& event, double trackFitLength,
                                  double trackFitDistanceFromStart) {
//***************************************************************
  constexpr Float_t kUnassigned = -999.0f;
  if (!neutral) return;
  neutral->AlignmentPandora = kUnassigned;
  neutral->AlignmentFit = kUnassigned;

  AnaAnnihilationVertexPD* vertex = neutral->AnnihilationVertex;
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* daughter1 = vertex->Particles[0];
  AnaParticlePD* daughter2 = vertex->Particles[1];
  if (!daughter1 || !daughter2) return;

  Float_t p1Mag = daughter1->Momentum;
  Float_t p2Mag = daughter2->Momentum;
  if (!(std::isfinite(p1Mag) && std::isfinite(p2Mag) && p1Mag > 0.0f && p2Mag > 0.0f)) {
    p1Mag = 1.0f;
    p2Mag = 1.0f;
  }

  TVector3 dirPandora1, dirPandora2, dirFit1, dirFit2;
  DaughterPandoraAndFitDirs(vertex, trackFitLength, trackFitDistanceFromStart, dirPandora1, dirPandora2, dirFit1,
                            dirFit2);

  if (dirPandora1.Mag2() <= 0.0 || dirPandora2.Mag2() <= 0.0) return;
  if (dirFit1.Mag2() <= 0.0 || dirFit2.Mag2() <= 0.0) return;

  const TVector3 pTotPandora = p1Mag * dirPandora1.Unit() + p2Mag * dirPandora2.Unit();
  const TVector3 pTotFit = p1Mag * dirFit1.Unit() + p2Mag * dirFit2.Unit();
  if (pTotPandora.Mag2() <= 0.0 || pTotFit.Mag2() <= 0.0) return;

  const Float_t* s = neutral->PositionStart;
  if (s[0] <= -900.0f || s[1] <= -900.0f || s[2] <= -900.0f) return;

  if (vertex->PositionPandora[0] > -900.0f && vertex->PositionPandora[1] > -900.0f &&
      vertex->PositionPandora[2] > -900.0f) {
    TVector3 uP(vertex->PositionPandora[0] - s[0], vertex->PositionPandora[1] - s[1],
                  vertex->PositionPandora[2] - s[2]);
    if (uP.Mag2() > 0.0) {
      neutral->AlignmentPandora = AngleBetweenUnitVectors(uP, pTotPandora);
    }
  }

  if (vertex->PositionFit[0] > -900.0f && vertex->PositionFit[1] > -900.0f && vertex->PositionFit[2] > -900.0f) {
    TVector3 uF(vertex->PositionFit[0] - s[0], vertex->PositionFit[1] - s[1], vertex->PositionFit[2] - s[2]);
    if (uF.Mag2() > 0.0) {
      neutral->AlignmentFit = AngleBetweenUnitVectors(uF, pTotFit);
    }
  }
}

//***************************************************************
void FillPositionPandora(AnaAnnihilationVertexPD* vertex) {
//***************************************************************
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* p1 = vertex->Particles[0];
  AnaParticlePD* p2 = vertex->Particles[1];
  if (!p1 || !p2) return;

  TVector3 x1(p1->PositionStart[0], p1->PositionStart[1], p1->PositionStart[2]);
  TVector3 x2(p2->PositionStart[0], p2->PositionStart[1], p2->PositionStart[2]);
  TVector3 d1(p1->DirectionStart[0], p1->DirectionStart[1], p1->DirectionStart[2]);
  TVector3 d2(p2->DirectionStart[0], p2->DirectionStart[1], p2->DirectionStart[2]);

  const double eps = 1e-10;
  if (d1.Mag2() < eps || d2.Mag2() < eps) {
    TVector3 fallback = 0.5 * (x1 + x2);
    vertex->PositionPandora[0] = fallback.X();
    vertex->PositionPandora[1] = fallback.Y();
    vertex->PositionPandora[2] = fallback.Z();
    vertex->ClosestPointPandora1[0] = x1.X();
    vertex->ClosestPointPandora1[1] = x1.Y();
    vertex->ClosestPointPandora1[2] = x1.Z();
    vertex->ClosestPointPandora2[0] = x2.X();
    vertex->ClosestPointPandora2[1] = x2.Y();
    vertex->ClosestPointPandora2[2] = x2.Z();
    vertex->MinimumDistancePandora = static_cast<float>((x1 - x2).Mag());
    return;
  }

  d1 = d1.Unit();
  d2 = d2.Unit();
  TVector3 w0 = x1 - x2;

  const double a = d1.Dot(d1);
  const double b = d1.Dot(d2);
  const double c = d2.Dot(d2);
  const double d = d1.Dot(w0);
  const double e = d2.Dot(w0);
  const double den = a * c - b * b;

  if (den < eps && den > -eps) {
    TVector3 fallback = 0.5 * (x1 + x2);
    vertex->PositionPandora[0] = fallback.X();
    vertex->PositionPandora[1] = fallback.Y();
    vertex->PositionPandora[2] = fallback.Z();
    vertex->ClosestPointPandora1[0] = x1.X();
    vertex->ClosestPointPandora1[1] = x1.Y();
    vertex->ClosestPointPandora1[2] = x1.Z();
    vertex->ClosestPointPandora2[0] = x2.X();
    vertex->ClosestPointPandora2[1] = x2.Y();
    vertex->ClosestPointPandora2[2] = x2.Z();
    vertex->MinimumDistancePandora = static_cast<float>((x1 - x2).Mag());
    return;
  }

  const double t1 = (b * e - c * d) / den;
  const double t2 = (a * e - b * d) / den;
  TVector3 c1 = x1 + t1 * d1;
  TVector3 c2 = x2 + t2 * d2;
  TVector3 intersection = 0.5 * (c1 + c2);
  const float minDistance = static_cast<float>((c1 - c2).Mag());

  vertex->PositionPandora[0] = intersection.X();
  vertex->PositionPandora[1] = intersection.Y();
  vertex->PositionPandora[2] = intersection.Z();
  vertex->ClosestPointPandora1[0] = c1.X();
  vertex->ClosestPointPandora1[1] = c1.Y();
  vertex->ClosestPointPandora1[2] = c1.Z();
  vertex->ClosestPointPandora2[0] = c2.X();
  vertex->ClosestPointPandora2[1] = c2.Y();
  vertex->ClosestPointPandora2[2] = c2.Z();
  vertex->MinimumDistancePandora = minDistance;
}

//***************************************************************
void FillPositionFit(AnaAnnihilationVertexPD* vertex, double trackFitLength, double trackFitDistanceFromStart) {
//***************************************************************
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* p1 = vertex->Particles[0];
  AnaParticlePD* p2 = vertex->Particles[1];
  if (!p1 || !p2) return;

  std::vector<double> fit1;
  std::vector<double> fit2;
  pdAnaUtils::ExtrapolateTrack(p1, fit1, trackFitLength, true, trackFitDistanceFromStart);
  pdAnaUtils::ExtrapolateTrack(p2, fit2, trackFitLength, true, trackFitDistanceFromStart);

  TVector3 s1(p1->PositionStart[0], p1->PositionStart[1], p1->PositionStart[2]);
  TVector3 s2(p2->PositionStart[0], p2->PositionStart[1], p2->PositionStart[2]);
  const bool fit1Valid = (fit1.size() >= 6 && fit1[0] > -900.0 && fit1[1] > -900.0 && fit1[2] > -900.0);
  const bool fit2Valid = (fit2.size() >= 6 && fit2[0] > -900.0 && fit2[1] > -900.0 && fit2[2] > -900.0);

  if (fit1Valid) {
    s1.SetXYZ(fit1[0], fit1[1], fit1[2]);
  }
  if (fit2Valid) {
    s2.SetXYZ(fit2[0], fit2[1], fit2[2]);
  }

  TVector3 d1(p1->DirectionStart[0], p1->DirectionStart[1], p1->DirectionStart[2]);
  TVector3 d2(p2->DirectionStart[0], p2->DirectionStart[1], p2->DirectionStart[2]);
  if (fit1Valid) d1.SetXYZ(fit1[3], fit1[4], fit1[5]);
  if (fit2Valid) d2.SetXYZ(fit2[3], fit2[4], fit2[5]);

  const double eps = 1e-10;
  if (d1.Mag2() < eps || d2.Mag2() < eps) {
    TVector3 fallback = 0.5 * (s1 + s2);
    vertex->PositionFit[0] = fallback.X();
    vertex->PositionFit[1] = fallback.Y();
    vertex->PositionFit[2] = fallback.Z();
    vertex->ClosestPointFit1[0] = s1.X();
    vertex->ClosestPointFit1[1] = s1.Y();
    vertex->ClosestPointFit1[2] = s1.Z();
    vertex->ClosestPointFit2[0] = s2.X();
    vertex->ClosestPointFit2[1] = s2.Y();
    vertex->ClosestPointFit2[2] = s2.Z();
    vertex->MinimumDistanceFit = static_cast<float>((s1 - s2).Mag());
    return;
  }
  d1 = d1.Unit();
  d2 = d2.Unit();

  std::vector<double> line1(6), line2(6);
  line1[0] = s1.X();     line1[1] = s1.Y();     line1[2] = s1.Z();
  line1[3] = d1.X();     line1[4] = d1.Y();     line1[5] = d1.Z();
  line2[0] = s2.X();     line2[1] = s2.Y();     line2[2] = s2.Z();
  line2[3] = d2.X();     line2[4] = d2.Y();     line2[5] = d2.Z();

  TVector3 closest1, closest2;
  const float minDistance = static_cast<float>(pdAnaUtils::FindClosestPointsBetweenLines(line1, line2, closest1, closest2));
  TVector3 fitPos = 0.5 * (closest1 + closest2);

  vertex->PositionFit[0] = fitPos.X();
  vertex->PositionFit[1] = fitPos.Y();
  vertex->PositionFit[2] = fitPos.Z();
  vertex->ClosestPointFit1[0] = closest1.X();
  vertex->ClosestPointFit1[1] = closest1.Y();
  vertex->ClosestPointFit1[2] = closest1.Z();
  vertex->ClosestPointFit2[0] = closest2.X();
  vertex->ClosestPointFit2[1] = closest2.Y();
  vertex->ClosestPointFit2[2] = closest2.Z();
  vertex->MinimumDistanceFit = minDistance;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> FilterVerticesByMinimumDistanceFit(std::vector<AnaAnnihilationVertexPD*>& vertices) {
//***************************************************************
  std::sort(vertices.begin(), vertices.end(),
            [](const AnaAnnihilationVertexPD* a, const AnaAnnihilationVertexPD* b) {
              return a->MinimumDistanceFit < b->MinimumDistanceFit;
            });

  std::unordered_set<AnaParticlePD*> usedParticles;
  std::vector<AnaAnnihilationVertexPD*> filtered;
  filtered.reserve(vertices.size());

  for (AnaAnnihilationVertexPD* vtx : vertices) {
    bool overlaps = false;
    for (AnaParticlePD* p : vtx->Particles) {
      if (usedParticles.find(p) != usedParticles.end()) {
        overlaps = true;
        break;
      }
    }

    if (overlaps) {
      delete vtx;
      continue;
    }

    filtered.push_back(vtx);
    for (AnaParticlePD* p : vtx->Particles) {
      usedParticles.insert(p);
    }
  }

  return filtered;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxDaughterDistance,
                                                           Int_t* nBeforeFiltering,
                                                           Int_t* nAfterFiltering) {
//***************************************************************
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart =
      ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
  const double annihilationDegeneracyRadius =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyRadius")
          : maxDaughterDistance;
  const double annihilationDegeneracyLineToVertexDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyLineToVertexDistance")
          : annihilationDegeneracyRadius;
  const double annihilationDegeneracyOriginSupportDistance =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDegeneracyOriginSupportDistance")
          : 0.5 * annihilationDegeneracyLineToVertexDistance;
  const int minCollectionHitsPerDaughter =
      ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexMinCollectionHits");
  const double annihilationProtonChi2NdfRejectBelow =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDaughterProtonChi2NdfRejectBelow")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDaughterProtonChi2NdfRejectBelow")
          : 100.0;
  const double annihilationPionChi2NdfRejectAbove =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDaughterPionChi2NdfRejectAbove")
          ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDaughterPionChi2NdfRejectAbove")
          : 25.0;
  const int annihilationUseEndpointCombinatorics =
      ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics")
          ? ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics")
          : 1;

  std::vector<AnaAnnihilationVertexPD*> reconstructedVertices;
  int vertexID = 0;

  for (int i = 0; i < nParts; ++i) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if (!daughter1) continue;
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) continue;
    if (daughter1->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;
    if (IsProtonLikeAndNotPionLike(daughter1, annihilationProtonChi2NdfRejectBelow,
                                   annihilationPionChi2NdfRejectAbove))
      continue;

    for (int j = i + 1; j < nParts; ++j) {
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if (!daughter2) continue;
      if (daughter1->ParentID != daughter2->ParentID) continue;
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) continue;
      if (daughter2->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;
      if (IsProtonLikeAndNotPionLike(daughter2, annihilationProtonChi2NdfRejectBelow,
                                     annihilationPionChi2NdfRejectAbove))
        continue;

      const AnnihilationEndpointPairing pairing =
          (annihilationUseEndpointCombinatorics != 0) ? BestAnnihilationPairEndpoint(daughter1, daughter2)
                                                      : StartStartAnnihilationPairing(daughter1, daughter2);
      if (!std::isfinite(pairing.distance) || pairing.distance > maxDaughterDistance) continue;

      AnaAnnihilationVertexPD* reconstructedVertex = new AnaAnnihilationVertexPD();
      reconstructedVertex->OriginalDistance = pairing.distance;
      reconstructedVertex->PairingEndpointCombo = pairing.combo;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);
      reconstructedVertex->NParticles = reconstructedVertex->Particles.size();

      // Pre-reversal fit for overlap filtering (MinimumDistanceFit); refined after ApplyAnnihilationPairingReversal.
      FillPositionPandora(reconstructedVertex);
      FillPositionFit(reconstructedVertex, trackFitLength, trackFitDistanceFromStart);
      reconstructedVertices.push_back(reconstructedVertex);
    }
  }

  if (nBeforeFiltering) {
    *nBeforeFiltering = static_cast<Int_t>(reconstructedVertices.size());
  }
  std::vector<AnaAnnihilationVertexPD*> filteredVertices = FilterVerticesByMinimumDistanceFit(reconstructedVertices);
  if (nAfterFiltering) {
    *nAfterFiltering = static_cast<Int_t>(filteredVertices.size());
  }

  for (AnaAnnihilationVertexPD* vertex : filteredVertices) {
    if (!vertex || vertex->Particles.size() < 2) continue;
    AnaParticlePD* d1 = vertex->Particles[0];
    AnaParticlePD* d2 = vertex->Particles[1];
    if (!d1 || !d2) continue;
    ApplyAnnihilationPairingReversal(vertex->PairingEndpointCombo, d1, d2);
    FillPositionPandora(vertex);
    FillPositionFit(vertex, trackFitLength, trackFitDistanceFromStart);
    vertex->Degeneracy =
        ComputeAnnihilationVertexDegeneracy(event, vertex, annihilationDegeneracyRadius,
                                           annihilationDegeneracyLineToVertexDistance,
                                           annihilationDegeneracyOriginSupportDistance, trackFitLength,
                                           trackFitDistanceFromStart);
  }

  for (AnaAnnihilationVertexPD* vertex : filteredVertices) {
    if (!vertex || vertex->Particles.size() < 2) continue;
    AnaParticlePD* daughter1 = vertex->Particles[0];
    AnaParticlePD* daughter2 = vertex->Particles[1];
    if (!daughter1 || !daughter2) continue;

    bool shouldAssignMomentum = true;
    if (ND::params().HasParameter("neutralKaonAnalysis.EnsureMomentumSignalOnly") &&
        ND::params().GetParameterI("neutralKaonAnalysis.EnsureMomentumSignalOnly") == 1) {
      const int signalCode = neutralKaonAnaUtils::GetSignalCategoryCodeForAnnihilationVertex(vertex, event);
      shouldAssignMomentum = (signalCode == 1 || signalCode == 5 || signalCode == 6);
    }

    DaughterMomentumDebugInfo daughter1Debug;
    DaughterMomentumDebugInfo daughter2Debug;
    vertex->JointK0sMomentumUsed = 0;
    vertex->JointK0sBestScore = -999.0f;
    vertex->JointK0sInvMassAtBest = -999.0f;
    vertex->JointK0sSigmaP1GeV = -999.0f;
    vertex->JointK0sSigmaP2GeV = -999.0f;
    vertex->JointK0sSigmaMEventGeV = -999.0f;
    vertex->JointK0sDmDp1 = -999.0f;
    vertex->JointK0sDmDp2 = -999.0f;
    vertex->JointK0sMomentumConstraintRatioR = -999.0f;
    vertex->JointK0sMomentumDedxChi2Degradation = -999.0f;
    vertex->JointK0sDebugClass = 0;

    if (shouldAssignMomentum) {
      const bool jointEnable = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumEnable") &&
                               ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMomentumEnable") == 1;

      vertex->Daughter1MomentumTLE = -999.0f;
      vertex->Daughter2MomentumTLE = -999.0f;
      vertex->Daughter1MomentumMCS = -999.0f;
      vertex->Daughter2MomentumMCS = -999.0f;
      vertex->Daughter1MomentumJoint = -999.0f;
      vertex->Daughter2MomentumJoint = -999.0f;
      bool jointApplied = false;
      if (jointEnable) {
        const double scanLmaxCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                                      ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                                      : 450.;
        const double scanStepCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepCm")
                                      ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepCm")
                                      : 1.;
        const int minInteriorHits =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
                ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
                : 15;
        const int skipFirst =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
                ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
                : 3;
        const int skipLast =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
                ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
                : 3;
        double dedxMinMeVcm =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
                : 0.5;
        double dedxMaxMeVcm =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
                : 5.0;
        const double pdfPathCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                                     ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                                     : 0.65;
        const double scanStepFineCm =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepFineCm")
                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepFineCm")
                : 0.;
        const double lowPMomentumRefineGeV =
            ND::params().HasParameter("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
                ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
                : 0.2;

        const double pMinGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                                   ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                                   : 0.05;
        const double pMaxGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                                   ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                                   : 2.0;
        const double pStepGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                                    ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                                    : 0.05;
        const double sigmaMassMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                          : 10.0;
        const double sigmaMassMinMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMinMeV")
                                           ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMinMeV")
                                           : 5.0;
        const double sigmaMassMaxMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMaxMeV")
                                           ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMaxMeV")
                                           : 50.0;
        const double penaltyScale = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                                        ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                                        : 1.0;
        const int refineFactor =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumRefineFactor")
                ? ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMomentumRefineFactor")
                : 2;
        const bool useMCS =
            !ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSEnable") ||
            ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMCSEnable") != 0;
        const double mcsWeight =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSWeight")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSWeight")
                : 1.0;
        pdJointK0sPionMomentum::MCSLikelihoodConfig mcsCfg;
        mcsCfg.radiationLengthCm =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSRadiationLengthCm")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSRadiationLengthCm")
                : 14.0;
        mcsCfg.minSegmentLengthCm =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSMinSegmentCm")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSMinSegmentCm")
                : 0.5;
        mcsCfg.theta0FloorRad =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSTheta0FloorRad")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSTheta0FloorRad")
                : 1e-6;
        mcsCfg.maxAbsDeltaThetaRad =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSMaxAbsDeltaThetaRad")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSMaxAbsDeltaThetaRad")
                : -1.0;
        mcsCfg.useDetectorSigma =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSUseDetectorSigma") &&
            ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMCSUseDetectorSigma") != 0;
        mcsCfg.detectorSigmaFloorRad =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaFloorRad")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaFloorRad")
                : 1e-6;
        mcsCfg.detectorSigmaA =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaA")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaA")
                : 0.0;
        mcsCfg.detectorSigmaC =
            ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaC")
                ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaC")
                : 0.0;
        constexpr double kK0sMassGeV = 0.497611;
        const double sigmaMassGeV = sigmaMassMeV * 1e-3;
        const double sigmaMassMinGeV = sigmaMassMinMeV * 1e-3;
        const double sigmaMassMaxGeV = sigmaMassMaxMeV * 1e-3;
        const double sigmaMassLoGeV = std::min(sigmaMassMinGeV, sigmaMassMaxGeV);
        const double sigmaMassHiGeV = std::max(sigmaMassMinGeV, sigmaMassMaxGeV);

        const bool useEventSigmaM =
            !ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") ||
            ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") != 0;

        std::vector<double> p1v, logL1, p2v, logL2;
        if (pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
                daughter1, scanLmaxCm, scanStepCm, minInteriorHits, skipFirst, skipLast, dedxMinMeVcm, dedxMaxMeVcm,
                pdfPathCm, p1v, logL1, scanStepFineCm, lowPMomentumRefineGeV) &&
            pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
                daughter2, scanLmaxCm, scanStepCm, minInteriorHits, skipFirst, skipLast, dedxMinMeVcm, dedxMaxMeVcm,
                pdfPathCm, p2v, logL2, scanStepFineCm, lowPMomentumRefineGeV)) {
          auto argmaxMomentumFromCurve = [](const std::vector<double>& pAxis, const std::vector<double>& logL) -> double {
            if (pAxis.empty() || pAxis.size() != logL.size()) return std::numeric_limits<double>::quiet_NaN();
            size_t im = 0;
            for (size_t i = 1; i < logL.size(); ++i) {
              if (logL[i] > logL[im]) im = i;
            }
            return pAxis[im];
          };

          const double p1_tle = argmaxMomentumFromCurve(p1v, logL1);
          const double p2_tle = argmaxMomentumFromCurve(p2v, logL2);
          if (std::isfinite(p1_tle) && p1_tle > 0.) vertex->Daughter1MomentumTLE = static_cast<Float_t>(p1_tle);
          if (std::isfinite(p2_tle) && p2_tle > 0.) vertex->Daughter2MomentumTLE = static_cast<Float_t>(p2_tle);

          std::vector<double> logL1Joint = logL1;
          std::vector<double> logL2Joint = logL2;
          if (useMCS) {
            std::vector<double> logLMcs1;
            std::vector<double> logLMcs2;
            const bool okMcs1 =
                pdJointK0sPionMomentum::BuildPionMCSLogLikelihoodVsMomentumCurve(*daughter1, p1v, mcsCfg, logLMcs1);
            const bool okMcs2 =
                pdJointK0sPionMomentum::BuildPionMCSLogLikelihoodVsMomentumCurve(*daughter2, p2v, mcsCfg, logLMcs2);
            if (okMcs1 && okMcs2 && logLMcs1.size() == logL1Joint.size() && logLMcs2.size() == logL2Joint.size()) {
              const double p1_mcs = argmaxMomentumFromCurve(p1v, logLMcs1);
              const double p2_mcs = argmaxMomentumFromCurve(p2v, logLMcs2);
              if (std::isfinite(p1_mcs) && p1_mcs > 0.) vertex->Daughter1MomentumMCS = static_cast<Float_t>(p1_mcs);
              if (std::isfinite(p2_mcs) && p2_mcs > 0.) vertex->Daughter2MomentumMCS = static_cast<Float_t>(p2_mcs);
              if (std::isfinite(mcsWeight) && mcsWeight > 0.) {
                for (size_t i = 0; i < logL1Joint.size(); ++i) logL1Joint[i] += mcsWeight * logLMcs1[i];
                for (size_t i = 0; i < logL2Joint.size(); ++i) logL2Joint[i] += mcsWeight * logLMcs2[i];
              }
            }
          }

          TVector3 dirP1, dirP2, dirFit1, dirFit2;
          DaughterPandoraAndFitDirs(vertex, trackFitLength, trackFitDistanceFromStart, dirP1, dirP2, dirFit1, dirFit2);

          double sigma_m_for_grid = sigmaMassGeV;
          JointK0sSigmaMEventDiagnostics sigmaDiag;
          if (useEventSigmaM &&
              pdJointK0sPionMomentum::ComputeSigmaMEventGeV(p1v, logL1, p2v, logL2, dirFit1, dirFit2, sigmaMassGeV,
                                                            sigmaDiag, sigmaMassLoGeV, sigmaMassHiGeV) &&
              std::isfinite(sigmaDiag.sigma_m_event_gev) && sigmaDiag.sigma_m_event_gev > 0.) {
            sigma_m_for_grid = sigmaDiag.sigma_m_event_gev;
            vertex->JointK0sSigmaP1GeV = static_cast<Float_t>(sigmaDiag.sigma_p1_gev);
            vertex->JointK0sSigmaP2GeV = static_cast<Float_t>(sigmaDiag.sigma_p2_gev);
            vertex->JointK0sDmDp1 = static_cast<Float_t>(sigmaDiag.dm_dp1);
            vertex->JointK0sDmDp2 = static_cast<Float_t>(sigmaDiag.dm_dp2);
          }
          if (std::isfinite(sigma_m_for_grid) && sigma_m_for_grid > 0.) {
            sigma_m_for_grid = std::max(sigma_m_for_grid, sigmaMassLoGeV);
            sigma_m_for_grid = std::min(sigma_m_for_grid, sigmaMassHiGeV);
          }
          vertex->JointK0sSigmaMEventGeV = static_cast<Float_t>(sigma_m_for_grid);

          const JointK0sPionMomentumGridResult jr = pdJointK0sPionMomentum::FitJointMomentaOnGrid(
              p1v, logL1Joint, p2v, logL2Joint, dirFit1, dirFit2, pMinGeV, pMaxGeV, pStepGeV, kK0sMassGeV, sigma_m_for_grid,
              penaltyScale, refineFactor);
          if (jr.ok && std::isfinite(static_cast<double>(jr.p1)) && std::isfinite(static_cast<double>(jr.p2)) &&
              jr.p1 > 0.f && jr.p2 > 0.f) {
            daughter1->Momentum = jr.p1;
            daughter2->Momentum = jr.p2;
            vertex->Daughter1MomentumJoint = jr.p1;
            vertex->Daughter2MomentumJoint = jr.p2;
            jointApplied = true;
            vertex->JointK0sMomentumUsed = 1;
            vertex->JointK0sBestScore = jr.bestScore;
            vertex->JointK0sInvMassAtBest = jr.invMassAtBest;
            vertex->JointK0sMomentumConstraintRatioR = jr.constraintRatioR;
            vertex->JointK0sMomentumDedxChi2Degradation = -999.0f;
            vertex->JointK0sDebugClass = 0;
            {
              TProfile* jointPionTmpl = nullptr;
              TGraph* jointPionKe = nullptr;
              if (pdMomShared::LoadPionTemplates(jointPionTmpl, jointPionKe,
                                                  "pionRangeEnergyGraph_joint_postfit_dedx")) {
                std::vector<double> jk_d1x, jk_d1r, jk_d2x, jk_d2r;
                const bool jk_ok1 = JointK0sPrepPionInterior(daughter1, minInteriorHits, skipFirst, skipLast,
                                                            dedxMinMeVcm, dedxMaxMeVcm, jk_d1x, jk_d1r);
                const bool jk_ok2 = JointK0sPrepPionInterior(daughter2, minInteriorHits, skipFirst, skipLast,
                                                            dedxMinMeVcm, dedxMaxMeVcm, jk_d2x, jk_d2r);
                if (jk_ok1 && jk_ok2) {
                  size_t im1 = 0;
                  size_t im2 = 0;
                  for (size_t i = 1; i < logL1.size(); ++i) {
                    if (logL1[i] > logL1[im1]) im1 = i;
                  }
                  for (size_t i = 1; i < logL2.size(); ++i) {
                    if (logL2[i] > logL2[im2]) im2 = i;
                  }
                  const double p1_marg = p1v[im1];
                  const double p2_marg = p2v[im2];
                  double e1m = 0.;
                  double e2m = 0.;
                  double e1j = 0.;
                  double e2j = 0.;
                  const bool sol1m =
                      JointK0sBestExtensionForMomentum(daughter1, jk_d1r, p1_marg, scanLmaxCm, jointPionKe, e1m);
                  const bool sol2m =
                      JointK0sBestExtensionForMomentum(daughter2, jk_d2r, p2_marg, scanLmaxCm, jointPionKe, e2m);
                  const bool sol1j = JointK0sBestExtensionForMomentum(daughter1, jk_d1r, static_cast<double>(jr.p1),
                                                                      scanLmaxCm, jointPionKe, e1j);
                  const bool sol2j = JointK0sBestExtensionForMomentum(daughter2, jk_d2r, static_cast<double>(jr.p2),
                                                                      scanLmaxCm, jointPionKe, e2j);
                  const double c1b = pdMomReconstruction::CalculateExtensionChi2(jk_d1x, jk_d1r, e1m, jointPionTmpl);
                  const double c2b = pdMomReconstruction::CalculateExtensionChi2(jk_d2x, jk_d2r, e2m, jointPionTmpl);
                  const double c1a = pdMomReconstruction::CalculateExtensionChi2(jk_d1x, jk_d1r, e1j, jointPionTmpl);
                  const double c2a = pdMomReconstruction::CalculateExtensionChi2(jk_d2x, jk_d2r, e2j, jointPionTmpl);
                  const auto chi2FiniteOk = [](double x) { return std::isfinite(x) && x < 9000.; };
                  const bool okBefore =
                      chi2FiniteOk(c1b) && chi2FiniteOk(c2b) && sol1m && sol2m;
                  const bool okAfter = chi2FiniteOk(c1a) && chi2FiniteOk(c2a) && sol1j && sol2j;
                  const double sumB = c1b + c2b;
                  const double sumA = c1a + c2a;
                  if (okBefore && okAfter) {
                    vertex->JointK0sMomentumDedxChi2Degradation = static_cast<Float_t>(sumA - sumB);
                  } else if (okAfter) {
                    // Independent-TLE χ² baseline unavailable: store χ² at joint extensions only.
                    vertex->JointK0sMomentumDedxChi2Degradation = static_cast<Float_t>(sumA);
                  }
                  const Float_t Rv = jr.constraintRatioR;
                  if (okBefore && okAfter) {
                    const double del = sumA - sumB;
                    if (Rv < 0.3f && del <= 0.) vertex->JointK0sDebugClass = 1;
                    else if (Rv >= 1.f) vertex->JointK0sDebugClass = 2;
                    else if (del > 5.) vertex->JointK0sDebugClass = 3;
                  } else if (Rv >= 1.f) {
                    vertex->JointK0sDebugClass = 2;
                  }
                }
              }
            }
            daughter1Debug.hasPreexistingMomentum = -1;
            daughter2Debug.hasPreexistingMomentum = -1;
            daughter1Debug.extensionAttempted = 1;
            daughter2Debug.extensionAttempted = 1;
            daughter1Debug.extensionValid = 1;
            daughter2Debug.extensionValid = 1;
            daughter1Debug.extensionChi2Ndf = jr.bestScore;
            daughter2Debug.extensionChi2Ndf = jr.bestScore;
            daughter1Debug.extensionNValidHits = CountValidCollectionPlaneHits(daughter1);
            daughter2Debug.extensionNValidHits = CountValidCollectionPlaneHits(daughter2);
            daughter1Debug.extensionDedxBias = -999.0f;
            daughter2Debug.extensionDedxBias = -999.0f;
            daughter1Debug.extensionDedxSigma = -999.0f;
            daughter2Debug.extensionDedxSigma = -999.0f;
            daughter1Debug.extensionDedxFitOk = 1;
            daughter2Debug.extensionDedxFitOk = 1;
            vertex->Daughter1MomentumMethod = kMomMethodJointK0sGrid;
            vertex->Daughter2MomentumMethod = kMomMethodJointK0sGrid;
          }
        }
      }

      if (!jointApplied) {
        vertex->Daughter1MomentumMethod = AssignPionMomentumFromResidualRange(daughter1, &daughter1Debug);
        vertex->Daughter2MomentumMethod = AssignPionMomentumFromResidualRange(daughter2, &daughter2Debug);
      }
      vertex->Daughter1HasPreexistingMomentum = daughter1Debug.hasPreexistingMomentum;
      vertex->Daughter2HasPreexistingMomentum = daughter2Debug.hasPreexistingMomentum;
      vertex->Daughter1ExtensionAttempted = daughter1Debug.extensionAttempted;
      vertex->Daughter2ExtensionAttempted = daughter2Debug.extensionAttempted;
      vertex->Daughter1ExtensionValid = daughter1Debug.extensionValid;
      vertex->Daughter2ExtensionValid = daughter2Debug.extensionValid;
      vertex->Daughter1ExtensionChi2Ndf = daughter1Debug.extensionChi2Ndf;
      vertex->Daughter2ExtensionChi2Ndf = daughter2Debug.extensionChi2Ndf;
      vertex->Daughter1ExtensionNValidHits = daughter1Debug.extensionNValidHits;
      vertex->Daughter2ExtensionNValidHits = daughter2Debug.extensionNValidHits;
      vertex->Daughter1ExtensionDedxBias = daughter1Debug.extensionDedxBias;
      vertex->Daughter2ExtensionDedxBias = daughter2Debug.extensionDedxBias;
      vertex->Daughter1ExtensionDedxSigma = daughter1Debug.extensionDedxSigma;
      vertex->Daughter2ExtensionDedxSigma = daughter2Debug.extensionDedxSigma;
      vertex->Daughter1ExtensionDedxFitOk = daughter1Debug.extensionDedxFitOk;
      vertex->Daughter2ExtensionDedxFitOk = daughter2Debug.extensionDedxFitOk;
    } else {
      vertex->Daughter1MomentumMethod = kMomMethodUnassigned;
      vertex->Daughter2MomentumMethod = kMomMethodUnassigned;
      vertex->Daughter1HasPreexistingMomentum = -1;
      vertex->Daughter2HasPreexistingMomentum = -1;
      vertex->Daughter1ExtensionAttempted = 0;
      vertex->Daughter2ExtensionAttempted = 0;
      vertex->Daughter1ExtensionValid = 0;
      vertex->Daughter2ExtensionValid = 0;
      vertex->Daughter1ExtensionChi2Ndf = -999.0f;
      vertex->Daughter2ExtensionChi2Ndf = -999.0f;
      vertex->Daughter1ExtensionNValidHits = 0;
      vertex->Daughter2ExtensionNValidHits = 0;
      vertex->Daughter1ExtensionDedxBias = -999.0f;
      vertex->Daughter2ExtensionDedxBias = -999.0f;
      vertex->Daughter1ExtensionDedxSigma = -999.0f;
      vertex->Daughter2ExtensionDedxSigma = -999.0f;
      vertex->Daughter1ExtensionDedxFitOk = -1;
      vertex->Daughter2ExtensionDedxFitOk = -1;
    }
    FillVertexKinematicsFromDaughters(vertex, trackFitLength, trackFitDistanceFromStart);
  }

  return filteredVertices;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance,
                                                     Int_t* nBeforeFiltering,
                                                     Int_t* nAfterFiltering) {
//***************************************************************
  return CreateVerticesCommon(event, maxDaughterDistance, nBeforeFiltering, nAfterFiltering);
}

} // namespace pdAnnihilationUtils