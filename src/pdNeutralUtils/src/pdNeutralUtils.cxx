#include "pdNeutralUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdCreationUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace pdNeutralUtils {
namespace {

/// Get the intermediate true particle between the reco parent and annihilation daughters.
/// The intermediate particle must be:
/// 1. The common true parent of the two annihilation vertex daughters
/// 2. A true daughter of the reco parent's true object
/// 3. Listed in its own parent's Daughters array
AnaTrueParticlePD* GetIntermediateTrueParticle(AnaEventB& event,
                                               AnaAnnihilationVertexPD* annihilationVtx,
                                               AnaParticlePD* recoParent) {
  if (!annihilationVtx || annihilationVtx->Particles.size() < 2) {
    return nullptr;
  }

  if (!recoParent || !recoParent->TrueObject) {
    return nullptr;
  }

  AnaParticlePD* daughter1 = annihilationVtx->Particles[0];
  AnaParticlePD* daughter2 = annihilationVtx->Particles[1];
  if (!daughter1 || !daughter2 || !daughter1->TrueObject || !daughter2->TrueObject) {
    return nullptr;
  }

  AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
  if (!trueDaughter1 || !trueDaughter2) {
    return nullptr;
  }

  // Verify both daughters have the same parent (the intermediate particle).
  if (trueDaughter1->ParentID <= 0 || trueDaughter1->ParentID != trueDaughter2->ParentID) {
    return nullptr;
  }

  AnaTrueParticlePD* intermediateParticle = pdAnaUtils::GetTrueParticle(&event, trueDaughter1->ParentID);
  if (!intermediateParticle || intermediateParticle->ID != trueDaughter1->ParentID) {
    return nullptr;
  }

  // Require exact daughter-parent consistency for the intermediate particle.
  bool hasDaughter1 = false;
  bool hasDaughter2 = false;
  for (size_t i = 0; i < intermediateParticle->Daughters.size(); ++i) {
    if (intermediateParticle->Daughters[i] == trueDaughter1->ID) hasDaughter1 = true;
    if (intermediateParticle->Daughters[i] == trueDaughter2->ID) hasDaughter2 = true;
  }
  if (!hasDaughter1 || !hasDaughter2) {
    return nullptr;
  }

  // Key validation: The intermediate particle must be a true daughter of the reco parent's true object.
  AnaTrueParticlePD* recoParentTrueObject = static_cast<AnaTrueParticlePD*>(recoParent->TrueObject);
  if (!recoParentTrueObject) {
    return nullptr;
  }

  if (intermediateParticle->ParentID != recoParentTrueObject->ID) {
    return nullptr;
  }

  // Verify that recoParent's true object lists the intermediate in its Daughters array.
  bool foundIntermediate = false;
  for (size_t i = 0; i < recoParentTrueObject->Daughters.size(); ++i) {
    if (recoParentTrueObject->Daughters[i] == intermediateParticle->ID) {
      foundIntermediate = true;
      break;
    }
  }
  if (!foundIntermediate) {
    return nullptr;
  }

  return intermediateParticle;
}

bool HasValidRecoPoint3(const Float_t* point) {
  if (!point) return false;
  return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]) &&
         point[0] > -900.f && point[1] > -900.f && point[2] > -900.f;
}

bool BuildRecoFitLine(const AnaParticlePD* particle,
                      double trackFitLength,
                      double trackFitDistanceFromStart,
                      std::vector<double>& fitLine) {
  fitLine.clear();
  if (!particle) return false;
  pdAnaUtils::ExtrapolateTrack(const_cast<AnaParticlePD*>(particle), fitLine, trackFitLength, true,
                               trackFitDistanceFromStart);
  const bool valid = (fitLine.size() >= 6 && std::isfinite(fitLine[0]) && std::isfinite(fitLine[1]) &&
                      std::isfinite(fitLine[2]) && std::isfinite(fitLine[3]) && std::isfinite(fitLine[4]) &&
                      std::isfinite(fitLine[5]) && fitLine[0] > -900.0 && fitLine[1] > -900.0 &&
                      fitLine[2] > -900.0);
  return valid;
}

bool IsInsideCreationCylinder(const TVector3& point,
                              const std::vector<double>& lineParams,
                              const TVector3& creationPos,
                              const TVector3& annihilationPos,
                              double cylinderRadius) {
  if (lineParams.size() < 6 || cylinderRadius <= 0.0) return false;
  TVector3 linePoint(lineParams[0], lineParams[1], lineParams[2]);
  TVector3 lineDir(lineParams[3], lineParams[4], lineParams[5]);
  if (lineDir.Mag2() <= 1e-10) return false;
  lineDir = lineDir.Unit();

  const double tPoint = (point - linePoint).Dot(lineDir);
  const double tCreation = (creationPos - linePoint).Dot(lineDir);
  const double tAnnihilation = (annihilationPos - linePoint).Dot(lineDir);
  const double tMin = std::min(tCreation, tAnnihilation);
  const double tMax = std::max(tCreation, tAnnihilation);
  if (tPoint < tMin || tPoint > tMax) return false;

  const TVector3 projected = linePoint + tPoint * lineDir;
  return (point - projected).Mag() <= cylinderRadius;
}

double EffectiveNeutralLengthCm(const AnaNeutralParticlePD* np) {
  if (!np) return 1e12;
  auto ok = [](Float_t x) { return x > 0.f && std::isfinite(x); };
  if (ok(np->LengthFit)) return static_cast<double>(np->LengthFit);
  if (ok(np->LengthPandora)) return static_cast<double>(np->LengthPandora);
  if (ok(np->Length)) return static_cast<double>(np->Length);
  return 1e12;
}

double NeutralBeamZAbsCosTheta(const AnaNeutralParticlePD* np) {
  if (!np) return -1.0;
  const float dz = np->DirectionStart[2];
  if (dz <= -900.f || !std::isfinite(dz)) return -1.0;
  return static_cast<double>(std::fabs(static_cast<double>(dz)));
}

int ParentRecoHitSum(const AnaParticlePD* p) {
  if (!p) return 0;
  return p->NHitsPerPlane[0] + p->NHitsPerPlane[1] + p->NHitsPerPlane[2];
}

/// sortedIdx: candidate indices best-to-worst; keys[i] is metric for candidate i.
/// tierPoints[k] is awarded to every candidate in the k-th distinct tier; tied values share a tier
/// and the next tier skips (e.g. two tied for 1st both get tierPoints[0], next group gets tierPoints[2] if 3 tiers).
void AddTieredLeaderboardPoints(std::vector<int>& scores,
                                const std::vector<size_t>& sortedIdx,
                                const std::vector<double>& keys,
                                const int* tierPoints,
                                int nTiers,
                                double eps) {
  const size_t n = sortedIdx.size();
  size_t pos = 0;
  for (int t = 0; t < nTiers && pos < n; ++t) {
    const double groupKey = keys[sortedIdx[pos]];
    size_t end = pos + 1;
    while (end < n && std::fabs(keys[sortedIdx[end]] - groupKey) <= eps) {
      ++end;
    }
    const int pts = tierPoints[t];
    for (size_t k = pos; k < end; ++k) {
      scores[sortedIdx[k]] += pts;
    }
    pos = end;
  }
}

/// When several parents share one annihilation vertex, keep one neutral: score from three matches
/// (length 5/3/1, alignment 3/2/1, parent hits 1 for top tier only). Ties share tier points and skip
/// intermediate tiers. Final tie-break: shortest effective length, then parent UniqueID.
size_t SelectWinnerNeutralCandidateIndex(const std::vector<AnaNeutralParticlePD*>& cands) {
  const size_t n = cands.size();
  if (n == 0) return 0;
  if (n == 1) return 0;

  std::vector<int> scores(n, 0);
  constexpr double kLenEps = 1e-5;
  constexpr double kBeamEps = 1e-9;

  std::vector<double> lenKey(n);
  std::vector<double> beamKey(n);
  std::vector<double> hitKey(n);
  for (size_t i = 0; i < n; ++i) {
    lenKey[i] = EffectiveNeutralLengthCm(cands[i]);
    beamKey[i] = NeutralBeamZAbsCosTheta(cands[i]);
    hitKey[i] = static_cast<double>(ParentRecoHitSum(cands[i] ? cands[i]->Parent : nullptr));
  }

  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return lenKey[a] < lenKey[b]; });
  {
    const int pts[] = {5, 3, 1};
    AddTieredLeaderboardPoints(scores, idx, lenKey, pts, 3, kLenEps);
  }

  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return beamKey[a] > beamKey[b]; });
  {
    const int pts[] = {3, 2, 1};
    AddTieredLeaderboardPoints(scores, idx, beamKey, pts, 3, kBeamEps);
  }

  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return hitKey[a] > hitKey[b]; });
  {
    const int pts[] = {1};
    AddTieredLeaderboardPoints(scores, idx, hitKey, pts, 1, 0.0);
  }

  size_t bestIdx = 0;
  int bestScore = scores[0];
  double bestLen = EffectiveNeutralLengthCm(cands[0]);
  int bestParentUid = (cands[0] && cands[0]->Parent) ? cands[0]->Parent->UniqueID : std::numeric_limits<int>::max();

  for (size_t i = 1; i < n; ++i) {
    const int sc = scores[i];
    const double len = EffectiveNeutralLengthCm(cands[i]);
    const int puid = (cands[i] && cands[i]->Parent) ? cands[i]->Parent->UniqueID : std::numeric_limits<int>::max();
    if (sc > bestScore) {
      bestIdx = i;
      bestScore = sc;
      bestLen = len;
      bestParentUid = puid;
    } else if (sc == bestScore) {
      if (len < bestLen - kLenEps) {
        bestIdx = i;
        bestLen = len;
        bestParentUid = puid;
      } else if (std::fabs(len - bestLen) <= kLenEps && puid < bestParentUid) {
        bestIdx = i;
        bestParentUid = puid;
      }
    }
  }
  return bestIdx;
}

} // namespace

//***************************************************************
std::vector<AnaNeutralParticlePD*> CreateNeutralsFromAnnihilationVertices(
    AnaEventB& event,
    const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices) {
//***************************************************************

  std::vector<AnaNeutralParticlePD*> neutralParticles;
  int neutralParticleID = 0;
  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart =
      ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
  const double creationVertexRadius =
      ND::params().HasParameter("neutralKaonAnalysis.CreationVertexRadius")
          ? ND::params().GetParameterD("neutralKaonAnalysis.CreationVertexRadius")
          : 5.0;
  const double creationVertexSecondParticleMaxProtonChi2Ndf =
      ND::params().HasParameter("neutralKaonAnalysis.CreationVertexSecondParticleMaxProtonChi2Ndf")
          ? ND::params().GetParameterD("neutralKaonAnalysis.CreationVertexSecondParticleMaxProtonChi2Ndf")
          : 100.0;
  const bool selectSingleNeutralPerAnnihilationVertex =
      !ND::params().HasParameter("neutralKaonAnalysis.SelectSingleNeutralPerAnnihilationVertex") ||
      ND::params().GetParameterD("neutralKaonAnalysis.SelectSingleNeutralPerAnnihilationVertex") != 0.0;

  for (AnaAnnihilationVertexPD* annihilationVtx : annihilationVertices) {
    if (!annihilationVtx) continue;

    AnaParticlePD* vertexDaughter1 = nullptr;
    AnaParticlePD* vertexDaughter2 = nullptr;
    if (annihilationVtx->Particles.size() >= 2) {
      vertexDaughter1 = annihilationVtx->Particles[0];
      vertexDaughter2 = annihilationVtx->Particles[1];
    }

    std::vector<AnaParticlePD*> parentParticles;
    TVector3 annihilationPos(-999.0, -999.0, -999.0);
    if (annihilationVtx->PositionFit[0] > -900.f &&
        annihilationVtx->PositionFit[1] > -900.f &&
        annihilationVtx->PositionFit[2] > -900.f) {
      annihilationPos.SetXYZ(annihilationVtx->PositionFit[0], annihilationVtx->PositionFit[1], annihilationVtx->PositionFit[2]);
    } else if (annihilationVtx->PositionPandora[0] > -900.f &&
               annihilationVtx->PositionPandora[1] > -900.f &&
               annihilationVtx->PositionPandora[2] > -900.f) {
      annihilationPos.SetXYZ(annihilationVtx->PositionPandora[0], annihilationVtx->PositionPandora[1], annihilationVtx->PositionPandora[2]);
    }

    if (annihilationPos.X() > -900.0 && annihilationPos.Y() > -900.0 && annihilationPos.Z() > -900.0) {
      for (Int_t i = 0; i < event.nParticles; ++i) {
        AnaParticlePD* candidateParent = static_cast<AnaParticlePD*>(event.Particles[i]);
        if (!candidateParent) continue;
        if (candidateParent == vertexDaughter1 || candidateParent == vertexDaughter2) continue;

        const bool hasValidEnd = candidateParent->PositionEnd[0] > -900.f &&
                                 candidateParent->PositionEnd[1] > -900.f &&
                                 candidateParent->PositionEnd[2] > -900.f;
        if (!hasValidEnd) continue;

        const TVector3 candidateEnd(candidateParent->PositionEnd[0],
                                    candidateParent->PositionEnd[1],
                                    candidateParent->PositionEnd[2]);
        if (candidateEnd.Z() >= annihilationPos.Z()) continue;

        parentParticles.push_back(candidateParent);
      }
    }

    if (parentParticles.empty()) continue;

    std::vector<AnaNeutralParticlePD*> vertexCandidates;
    vertexCandidates.reserve(parentParticles.size());

    for (size_t ip = 0; ip < parentParticles.size(); ++ip) {
      AnaParticlePD* parentParticle = parentParticles[ip];
      if (!parentParticle) continue;

      const bool hasValidRecoParentEnd = parentParticle->PositionEnd[0] > -900.f &&
                                         parentParticle->PositionEnd[1] > -900.f &&
                                         parentParticle->PositionEnd[2] > -900.f;
      if (!hasValidRecoParentEnd) continue;

      AnaCreationVertexPD* creationVtx = new AnaCreationVertexPD();
      creationVtx->BeamParticle = parentParticle;
      creationVtx->SecondParticle = nullptr;

      TVector3 correctedBeamEnd(-999.0, -999.0, -999.0);
      bool hasProjectedBeamEnd = false;
      if (parentParticle) {
        double fitDistanceFromEndCm = 10.0;
        if (ND::params().HasParameter("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd")) {
          fitDistanceFromEndCm = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd");
        }
        hasProjectedBeamEnd = pdCreationUtils::ProjectBeamTailOntoStartDirection(parentParticle,
                                                                                  fitDistanceFromEndCm,
                                                                                  correctedBeamEnd);
      }

      const bool hasValidProjectedBeamEnd = hasProjectedBeamEnd && correctedBeamEnd.X() > -900.0 &&
                                            correctedBeamEnd.Y() > -900.0 && correctedBeamEnd.Z() > -900.0;

      const Float_t parentEndX = (hasValidProjectedBeamEnd ? static_cast<Float_t>(correctedBeamEnd.X())
                                                            : (parentParticle ? parentParticle->PositionEnd[0] : -999.f));
      const Float_t parentEndY = (hasValidProjectedBeamEnd ? static_cast<Float_t>(correctedBeamEnd.Y())
                                                            : (parentParticle ? parentParticle->PositionEnd[1] : -999.f));
      const Float_t parentEndZ = (hasValidProjectedBeamEnd ? static_cast<Float_t>(correctedBeamEnd.Z())
                                                            : (parentParticle ? parentParticle->PositionEnd[2] : -999.f));

      const bool hasValidAnnihilationPandoraPos = annihilationVtx->PositionPandora[0] > -900.f &&
                                                  annihilationVtx->PositionPandora[1] > -900.f &&
                                                  annihilationVtx->PositionPandora[2] > -900.f;
      if (hasValidAnnihilationPandoraPos && parentEndZ >= annihilationVtx->PositionPandora[2]) {
        continue;
      }

      if (parentParticle) {
        creationVtx->UniqueID = parentParticle->UniqueID;
        creationVtx->Particles.push_back(parentParticle);
        creationVtx->NParticles = static_cast<Int_t>(creationVtx->Particles.size());
        creationVtx->PositionPandora[0] = parentEndX;
        creationVtx->PositionPandora[1] = parentEndY;
        creationVtx->PositionPandora[2] = parentEndZ;
        creationVtx->Position[0] = parentEndX;
        creationVtx->Position[1] = parentEndY;
        creationVtx->Position[2] = parentEndZ;
      }

      const bool hasValidCreationPos = parentEndX > -900.f && parentEndY > -900.f && parentEndZ > -900.f;
      const bool hasValidAnnihilationPosForCylinder = annihilationPos.X() > -900.0 && annihilationPos.Y() > -900.0 &&
                                                      annihilationPos.Z() > -900.0;
      if (parentParticle && hasValidCreationPos && hasValidAnnihilationPosForCylinder && creationVertexRadius > 0.0 &&
          trackFitLength > 0.0 && trackFitDistanceFromStart >= 0.0) {
        std::vector<double> beamLineParams;
        if (BuildRecoFitLine(parentParticle, trackFitLength, trackFitDistanceFromStart, beamLineParams)) {
          const TVector3 creationPos(parentEndX, parentEndY, parentEndZ);
          Float_t bestMinDistance = 1e9f;
          std::vector<double> bestSecondLineParams;
          TVector3 bestClosestBeam(-999.0, -999.0, -999.0);
          TVector3 bestClosestSecond(-999.0, -999.0, -999.0);
          AnaParticlePD* bestSecondParticle = nullptr;
          Float_t bestProtonScore = -999.f;
          Float_t bestDistanceScore = -999.f;

          for (Int_t ipart = 0; ipart < event.nParticles; ++ipart) {
            AnaParticlePD* candidateSecond = static_cast<AnaParticlePD*>(event.Particles[ipart]);
            if (!candidateSecond) continue;
            if (candidateSecond == parentParticle || candidateSecond == vertexDaughter1 || candidateSecond == vertexDaughter2) continue;
            if (!HasValidRecoPoint3(candidateSecond->PositionStart)) continue;

            const Float_t protonScore = pdAnaUtils::Chi2PIDChi2PerHit(candidateSecond, 2212);
            if (!(protonScore > -900.f)) continue;
            if (!(std::isfinite(protonScore) && protonScore <= creationVertexSecondParticleMaxProtonChi2Ndf)) continue;

            const TVector3 candidateStart(candidateSecond->PositionStart[0],
                                          candidateSecond->PositionStart[1],
                                          candidateSecond->PositionStart[2]);
            if (!IsInsideCreationCylinder(candidateStart, beamLineParams, creationPos, annihilationPos, creationVertexRadius)) {
              continue;
            }

            std::vector<double> secondLineParams;
            if (!BuildRecoFitLine(candidateSecond, trackFitLength, trackFitDistanceFromStart, secondLineParams)) continue;

            TVector3 closestBeam;
            TVector3 closestSecond;
            const Float_t minDistance =
                static_cast<Float_t>(pdAnaUtils::FindClosestPointsBetweenLines(beamLineParams, secondLineParams,
                                                                               closestBeam, closestSecond));
            if (!(std::isfinite(minDistance) && minDistance >= 0.f)) continue;
            if (minDistance >= bestMinDistance) continue;

            bestMinDistance = minDistance;
            bestSecondParticle = candidateSecond;
            bestSecondLineParams = secondLineParams;
            bestClosestBeam = closestBeam;
            bestClosestSecond = closestSecond;
            bestProtonScore = protonScore;
            bestDistanceScore = static_cast<Float_t>((candidateStart - creationPos).Mag());
          }

          if (bestSecondParticle) {
            creationVtx->SecondParticle = bestSecondParticle;
            creationVtx->ProtonScore = bestProtonScore;
            creationVtx->DistanceScore = bestDistanceScore;
            creationVtx->MinDistanceScore = bestMinDistance;
            creationVtx->ClosestPointBeam[0] = static_cast<Float_t>(bestClosestBeam.X());
            creationVtx->ClosestPointBeam[1] = static_cast<Float_t>(bestClosestBeam.Y());
            creationVtx->ClosestPointBeam[2] = static_cast<Float_t>(bestClosestBeam.Z());
            creationVtx->ClosestPointSecond[0] = static_cast<Float_t>(bestClosestSecond.X());
            creationVtx->ClosestPointSecond[1] = static_cast<Float_t>(bestClosestSecond.Y());
            creationVtx->ClosestPointSecond[2] = static_cast<Float_t>(bestClosestSecond.Z());
            creationVtx->Position[0] = static_cast<Float_t>(0.5 * (bestClosestBeam.X() + bestClosestSecond.X()));
            creationVtx->Position[1] = static_cast<Float_t>(0.5 * (bestClosestBeam.Y() + bestClosestSecond.Y()));
            creationVtx->Position[2] = static_cast<Float_t>(0.5 * (bestClosestBeam.Z() + bestClosestSecond.Z()));
            creationVtx->PositionPandora[0] = creationVtx->Position[0];
            creationVtx->PositionPandora[1] = creationVtx->Position[1];
            creationVtx->PositionPandora[2] = creationVtx->Position[2];
            creationVtx->Particles.push_back(bestSecondParticle);
            creationVtx->NParticles = static_cast<Int_t>(creationVtx->Particles.size());
            creationVtx->FittedLineParams.clear();
            creationVtx->FittedLineParams.push_back(beamLineParams);
            creationVtx->FittedLineParams.push_back(bestSecondLineParams);

            // Only fill second-particle momentum for neutrals where annihilation-daughter
            // momentum assignment was enabled (same subset policy as annihilation daughters).
            if (annihilationVtx &&
                annihilationVtx->Daughter1MomentumMethod > -1 &&
                annihilationVtx->Daughter2MomentumMethod > -1) {
              pdAnnihilationUtils::AssignProtonMomentumFromResidualRange(bestSecondParticle);
            }
          }
        }
      }

      if (creationVtx) {
        creationVtx->Degeneracy = pdAnnihilationUtils::ComputeCreationVertexDegeneracy(
            event, creationVtx, annihilationVtx, -1);
      }

      AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();
      neutralParticle->AnnihilationVertex = annihilationVtx;
      neutralParticle->CreationVertex = creationVtx;
      neutralParticle->Parent = creationVtx->BeamParticle;

      const Float_t creationStartX = creationVtx->Position[0];
      const Float_t creationStartY = creationVtx->Position[1];
      const Float_t creationStartZ = creationVtx->Position[2];
      const bool hasValidParentEnd = parentParticle &&
                                     creationStartX > -900.f &&
                                     creationStartY > -900.f &&
                                     creationStartZ > -900.f;
      for (int i = 0; i < 3; ++i) {
        if (!hasValidParentEnd) {
          neutralParticle->PositionStart[i] = -999.f;
        } else if (i == 0) {
          neutralParticle->PositionStart[i] = creationStartX;
        } else if (i == 1) {
          neutralParticle->PositionStart[i] = creationStartY;
        } else {
          neutralParticle->PositionStart[i] = creationStartZ;
        }
        neutralParticle->PositionEnd[i] = annihilationVtx->PositionPandora[i];
        neutralParticle->DirectionEnd[i] = -999.f;
      }
      neutralParticle->PositionStart[3] = (parentParticle ? parentParticle->PositionEnd[3] : -999.f);
      neutralParticle->PositionEnd[3] = -999.f;

      const bool hasValidAnnihilationFitPos = annihilationVtx->PositionFit[0] > -900.f &&
                                              annihilationVtx->PositionFit[1] > -900.f &&
                                              annihilationVtx->PositionFit[2] > -900.f;

      neutralParticle->LengthPandora = -999.f;
      neutralParticle->LengthFit = -999.f;

      if (hasValidParentEnd && hasValidAnnihilationPandoraPos) {
        const TVector3 startPos(neutralParticle->PositionStart[0], neutralParticle->PositionStart[1], neutralParticle->PositionStart[2]);
        const TVector3 endPos(neutralParticle->PositionEnd[0], neutralParticle->PositionEnd[1], neutralParticle->PositionEnd[2]);
        TVector3 direction = endPos - startPos;
        neutralParticle->LengthPandora = static_cast<Float_t>(direction.Mag());
        neutralParticle->Length = neutralParticle->LengthPandora;
        if (direction.Mag2() > 0.) {
          direction = direction.Unit();
          neutralParticle->DirectionStart[0] = direction.X();
          neutralParticle->DirectionStart[1] = direction.Y();
          neutralParticle->DirectionStart[2] = direction.Z();
        } else {
          neutralParticle->DirectionStart[0] = -999.f;
          neutralParticle->DirectionStart[1] = -999.f;
          neutralParticle->DirectionStart[2] = -999.f;
        }
      } else {
        neutralParticle->Length = -999.f;
        neutralParticle->DirectionStart[0] = -999.f;
        neutralParticle->DirectionStart[1] = -999.f;
        neutralParticle->DirectionStart[2] = -999.f;
      }

      if (hasValidParentEnd && hasValidAnnihilationFitPos) {
        const TVector3 startPos(neutralParticle->PositionStart[0], neutralParticle->PositionStart[1], neutralParticle->PositionStart[2]);
        const TVector3 endPosFit(annihilationVtx->PositionFit[0], annihilationVtx->PositionFit[1], annihilationVtx->PositionFit[2]);
        neutralParticle->LengthFit = static_cast<Float_t>((endPosFit - startPos).Mag());
      }

      pdAnnihilationUtils::FillNeutralParticleAlignment(neutralParticle, event, trackFitLength,
                                                        trackFitDistanceFromStart);

      AnaTrueParticlePD* intermediateParticle = GetIntermediateTrueParticle(event, annihilationVtx, parentParticle);
      neutralParticle->TrueObject = intermediateParticle;

      vertexCandidates.push_back(neutralParticle);
    }

    if (vertexCandidates.empty()) continue;

    if (selectSingleNeutralPerAnnihilationVertex) {
      const size_t winnerIdx = (vertexCandidates.size() == 1)
                                   ? 0
                                   : SelectWinnerNeutralCandidateIndex(vertexCandidates);
      for (size_t ic = 0; ic < vertexCandidates.size(); ++ic) {
        if (ic == winnerIdx) {
          vertexCandidates[ic]->UniqueID = neutralParticleID++;
          neutralParticles.push_back(vertexCandidates[ic]);
        } else {
          delete vertexCandidates[ic]->CreationVertex;
          delete vertexCandidates[ic];
        }
      }
    } else {
      for (size_t ic = 0; ic < vertexCandidates.size(); ++ic) {
        vertexCandidates[ic]->UniqueID = neutralParticleID++;
        neutralParticles.push_back(vertexCandidates[ic]);
      }
    }
  }

  return neutralParticles;
}

} // namespace pdNeutralUtils