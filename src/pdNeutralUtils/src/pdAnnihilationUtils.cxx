#include "pdAnnihilationUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "neutralKaonAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace pdAnnihilationUtils {
namespace {

bool HasValidStartPosition(const AnaParticlePD* particle) {
  if (!particle) return false;
  return (particle->PositionStart[0] > -900.0 &&
          particle->PositionStart[1] > -900.0 &&
          particle->PositionStart[2] > -900.0);
}

bool HasValidPosition3(const Float_t pos[3]) {
  if (!pos) return false;
  return (pos[0] > -900.0f && pos[1] > -900.0f && pos[2] > -900.0f);
}

TVector3 GetAnnihilationVertexPositionForDegeneracy(const AnaAnnihilationVertexPD* vertex) {
  if (vertex && HasValidPosition3(vertex->PositionFit)) {
    return TVector3(vertex->PositionFit[0], vertex->PositionFit[1], vertex->PositionFit[2]);
  }
  if (vertex && HasValidPosition3(vertex->PositionPandora)) {
    return TVector3(vertex->PositionPandora[0], vertex->PositionPandora[1], vertex->PositionPandora[2]);
  }
  return TVector3(-999.0, -999.0, -999.0);
}

bool IsVertexDaughter(const AnaAnnihilationVertexPD* vertex, const AnaParticlePD* candidate) {
  if (!vertex || !candidate) return false;
  for (AnaParticlePD* daughter : vertex->Particles) {
    if (daughter == candidate) return true;
  }
  return false;
}

bool IsProtonLikeAndNotPionLike(const AnaParticlePD* particle) {
  if (!particle) return false;
  if (particle->Chi2ndf <= 0.f || particle->Chi2Proton <= 0.f) return false;

  const double protonChi2Ndf = particle->Chi2Proton / particle->Chi2ndf;
  if (protonChi2Ndf < 100.0) return true;

  const std::pair<double, int> pionPid = pdAnaUtils::Chi2PID(*particle, 211);
  if (pionPid.first < 0.0 || pionPid.second <= 0) return false;

  const double pionChi2Ndf = pionPid.first / pionPid.second;
  return (pionChi2Ndf > 25.0);
}

enum DaughterMomentumMethod {
  kMomMethodUnassigned = -1,
  kMomMethodPionRangeStopping = 0,
  kMomMethodCalorimetric = 1,
  kMomMethodFailed = 2,
  kMomMethodFreeRangeML = 3
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

Int_t AssignPionMomentumFromResidualRange(AnaParticlePD* particle, DaughterMomentumDebugInfo* debugInfo = nullptr) {
  if (!particle) return kMomMethodUnassigned;
  const Int_t nCollHits = static_cast<Int_t>(particle->Hits[2].size());
  const Int_t attempted = (nCollHits >= 4) ? 1 : 0;
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

  double truncMinRR = 0.;
  double truncFrac = 0.;
  if (ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxLandauTruncMinRRCm")) {
    truncMinRR = ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxLandauTruncMinRRCm");
  }
  if (ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxLandauTruncFraction")) {
    truncFrac = ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxLandauTruncFraction");
  }
  const pdAnaUtils::DEdxFreeRangeFitResult fit = pdAnaUtils::GetdEdxLikelihoodFreeRangeFit(
      particle, 211, 500., 0.5, truncMinRR, truncFrac);

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
                                          double annihilationVertexRadius) {
  if (!vertex || annihilationVertexRadius <= 0.0) return 0;

  const TVector3 vertexPos = GetAnnihilationVertexPositionForDegeneracy(vertex);
  if (vertexPos.X() < -900.0 || vertexPos.Y() < -900.0 || vertexPos.Z() < -900.0) return 0;

  Int_t degeneracy = 0;
  for (Int_t p = 0; p < event.nParticles; ++p) {
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[p]);
    if (!HasValidStartPosition(particle)) continue;
    if (IsVertexDaughter(vertex, particle)) continue;

    const TVector3 startPos(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
    if ((startPos - vertexPos).Mag() <= annihilationVertexRadius) {
      ++degeneracy;
    }
  }

  return degeneracy;
}

} // namespace

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
  const int minCollectionHitsPerDaughter =
      ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexMinCollectionHits");

  std::vector<AnaAnnihilationVertexPD*> reconstructedVertices;
  int vertexID = 0;

  for (int i = 0; i < nParts; ++i) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if (!daughter1) continue;
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) continue;
    if (daughter1->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;
    if (IsProtonLikeAndNotPionLike(daughter1)) continue;

    for (int j = i + 1; j < nParts; ++j) {
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if (!daughter2) continue;
      if (daughter1->ParentID != daughter2->ParentID) continue;
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) continue;
      if (daughter2->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;
      if (IsProtonLikeAndNotPionLike(daughter2)) continue;

      TVector3 s1(daughter1->PositionStart[0], daughter1->PositionStart[1], daughter1->PositionStart[2]);
      TVector3 s2(daughter2->PositionStart[0], daughter2->PositionStart[1], daughter2->PositionStart[2]);
      float distance = (s1 - s2).Mag();
      if (distance > maxDaughterDistance) continue;

      AnaAnnihilationVertexPD* reconstructedVertex = new AnaAnnihilationVertexPD();
      reconstructedVertex->OriginalDistance = distance;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);
      reconstructedVertex->NParticles = reconstructedVertex->Particles.size();

      FillPositionPandora(reconstructedVertex);
      FillPositionFit(reconstructedVertex, trackFitLength, trackFitDistanceFromStart);
      reconstructedVertex->Degeneracy =
          ComputeAnnihilationVertexDegeneracy(event, reconstructedVertex, maxDaughterDistance);
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
    if (shouldAssignMomentum) {
      vertex->Daughter1MomentumMethod = AssignPionMomentumFromResidualRange(daughter1, &daughter1Debug);
      vertex->Daughter2MomentumMethod = AssignPionMomentumFromResidualRange(daughter2, &daughter2Debug);
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