#include "pdNeutralUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdCreationUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <unordered_map>

namespace pdNeutralUtils {
namespace {

AnaTrueParticlePD* GetCommonTrueParent(AnaEventB& event, AnaAnnihilationVertexPD* annihilationVtx) {
  if (!annihilationVtx || annihilationVtx->Particles.size() < 2) {
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

  if (trueDaughter1->ParentID <= 0 || trueDaughter1->ParentID != trueDaughter2->ParentID) {
    return nullptr;
  }

  AnaTrueParticlePD* trueParent = pdAnaUtils::GetTrueParticle(&event, trueDaughter1->ParentID);
  if (!trueParent || trueParent->ID != trueDaughter1->ParentID) {
    return nullptr;
  }

  // Require exact daughter-parent consistency in truth.
  bool hasDaughter1 = false;
  bool hasDaughter2 = false;
  for (size_t i = 0; i < trueParent->Daughters.size(); ++i) {
    if (trueParent->Daughters[i] == trueDaughter1->ID) hasDaughter1 = true;
    if (trueParent->Daughters[i] == trueDaughter2->ID) hasDaughter2 = true;
  }
  if (!hasDaughter1 || !hasDaughter2) {
    return nullptr;
  }

  return trueParent;
}

} // namespace

//***************************************************************
std::vector<AnaNeutralParticlePD*> CreateNeutralsFromAnnihilationVertices(
    AnaEventB& event,
    const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices) {
//***************************************************************

  std::vector<AnaNeutralParticlePD*> neutralParticles;
  int neutralParticleID = 0;
  std::unordered_map<Int_t, AnaParticlePD*> particleByUniqueID;
  for (Int_t i = 0; i < event.nParticles; ++i) {
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[i]);
    if (particle) {
      particleByUniqueID[particle->UniqueID] = particle;
    }
  }

  for (AnaAnnihilationVertexPD* annihilationVtx : annihilationVertices) {
    if (!annihilationVtx) continue;

    AnaParticlePD* parentParticle = nullptr;
    if (annihilationVtx->Particles.size() >= 2) {
      AnaParticlePD* daughter1 = annihilationVtx->Particles[0];
      AnaParticlePD* daughter2 = annihilationVtx->Particles[1];
      if (daughter1 && daughter2 &&
          daughter1->ParentID > 0 &&
          daughter1->ParentID == daughter2->ParentID) {
        auto parentIt = particleByUniqueID.find(daughter1->ParentID);
        if (parentIt != particleByUniqueID.end()) {
          parentParticle = parentIt->second;
        }
      }
    }

    // Do not create a neutral candidate unless the reco parent exists and has a valid end position.
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

    const bool hasValidAnnihilationPandoraPos = annihilationVtx->PositionPandora[0] > -900.f &&
                                                annihilationVtx->PositionPandora[1] > -900.f &&
                                                annihilationVtx->PositionPandora[2] > -900.f;
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

    AnaTrueParticlePD* trueParentParticle = GetCommonTrueParent(event, annihilationVtx);
    neutralParticle->TrueObject = trueParentParticle;

    neutralParticles.push_back(neutralParticle);
  }

  return neutralParticles;
}

} // namespace pdNeutralUtils