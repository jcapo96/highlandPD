#include "pdNeutralUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdCreationUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"

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

} // namespace

//***************************************************************
std::vector<AnaNeutralParticlePD*> CreateNeutralsFromAnnihilationVertices(
    AnaEventB& event,
    const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices) {
//***************************************************************

  std::vector<AnaNeutralParticlePD*> neutralParticles;
  int neutralParticleID = 0;

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

      AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();
      neutralParticle->UniqueID = neutralParticleID++;
      neutralParticle->AnnihilationVertex = annihilationVtx;
      neutralParticle->CreationVertex = creationVtx;
      neutralParticle->Parent = creationVtx->BeamParticle;

      const bool hasValidParentEnd = parentParticle &&
                                     parentEndX > -900.f &&
                                     parentEndY > -900.f &&
                                     parentEndZ > -900.f;
      for (int i = 0; i < 3; ++i) {
        if (!hasValidParentEnd) {
          neutralParticle->PositionStart[i] = -999.f;
        } else if (i == 0) {
          neutralParticle->PositionStart[i] = parentEndX;
        } else if (i == 1) {
          neutralParticle->PositionStart[i] = parentEndY;
        } else {
          neutralParticle->PositionStart[i] = parentEndZ;
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

      const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
      const double trackFitDistanceFromStart =
          ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
      pdAnnihilationUtils::FillNeutralParticleAlignment(neutralParticle, event, trackFitLength,
                                                        trackFitDistanceFromStart);

      AnaTrueParticlePD* intermediateParticle = GetIntermediateTrueParticle(event, annihilationVtx, parentParticle);
      neutralParticle->TrueObject = intermediateParticle;

      neutralParticles.push_back(neutralParticle);
    }
  }

  return neutralParticles;
}

} // namespace pdNeutralUtils