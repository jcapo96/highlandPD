#include "pdNeutralUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdAnnihilationUtils.hxx"
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

    AnaCreationVertexPD* creationVtx = new AnaCreationVertexPD();
    creationVtx->BeamParticle = parentParticle;
    creationVtx->SecondParticle = nullptr;
    if (parentParticle) {
      creationVtx->UniqueID = parentParticle->UniqueID;
      creationVtx->Particles.push_back(parentParticle);
      creationVtx->NParticles = static_cast<Int_t>(creationVtx->Particles.size());
      creationVtx->PositionPandora[0] = parentParticle->PositionEnd[0];
      creationVtx->PositionPandora[1] = parentParticle->PositionEnd[1];
      creationVtx->PositionPandora[2] = parentParticle->PositionEnd[2];
      creationVtx->Position[0] = parentParticle->PositionEnd[0];
      creationVtx->Position[1] = parentParticle->PositionEnd[1];
      creationVtx->Position[2] = parentParticle->PositionEnd[2];
    }

    AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();
    neutralParticle->UniqueID = neutralParticleID++;
    neutralParticle->AnnihilationVertex = annihilationVtx;
    neutralParticle->CreationVertex = creationVtx;
    neutralParticle->Parent = creationVtx->BeamParticle;

    const bool hasValidParentEnd = parentParticle &&
                                   parentParticle->PositionEnd[0] > -900.f &&
                                   parentParticle->PositionEnd[1] > -900.f &&
                                   parentParticle->PositionEnd[2] > -900.f;
    for (int i = 0; i < 3; ++i) {
      neutralParticle->PositionStart[i] = hasValidParentEnd ? parentParticle->PositionEnd[i] : -999.f;
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

