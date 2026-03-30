#include "neutralKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdNeutralUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdBaseSelection.hxx"
#include "pdDataClasses.hxx"
#include "TVector3.h"
#include <algorithm>
#include <set>
#include <sstream>

//********************************************************************
neutralKaonSelection::neutralKaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void neutralKaonSelection::DefineSteps(){
//********************************************************************

  //basic cuts → see pdBaseAnalysis/src/pdBaseSelection
  AddStep(StepBase::kAction, "find Pandora track", new FindBeamTrackAction(), true);
  AddStep(StepBase::kCut   , "beam track in TPC"     , new BeamTrackExistsCut() , true);
  // AddStep(StepBase::kCut, "beam pdg filter", new BeamPDGCut(321), true);

  // Preliminary K0 Selection Steps - Look for neutral particle candidates
  AddStep(StepBase::kAction, "find all particles"   , new FindAllParticlesAction(), true);
  AddStep(StepBase::kCut   , "enough particles" , new HasEnoughParticlesCut(), true);
  AddStep(StepBase::kAction, "find neutral candidates", new FindNeutralCandidatesAction(), true);
  AddStep(StepBase::kCut   , "has neutral candidates" , new HasNeutralCandidatesCut(), true);

  // Minimal branch: keep reconstruction-stage candidate building only.
  // Neutral-particle purity cuts are intentionally disabled.

  // // Split the selection in branches, one for each possible candidate
  // AddSplit(neutralKaonAnalysisConstants::NMAXSAVEDCANDIDATES);

  // // Add branch-specific cuts for each candidate
  // for(int i = 0; i < (int)neutralKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
  //   AddStep(i, StepBase::kCut, "vertex daughters are daughters of vertex parent", new VtxDaughtersAreDaughtersOfVtxParentCut(), true);
  //   // AddStep(i, StepBase::kAction, "check K0 in truth"     , new CheckK0InTruthAction(), true);
  // }

  // // Set the branch aliases to the different branches
  // for(int i = 0; i < (int)neutralKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
  //   std::stringstream ssi;
  //   ssi << i;
  //   SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  // }
  SetBranchAlias(0,"trunk");
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
void neutralKaonSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB);

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void neutralKaonSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

//********************************************************************
bool FindAllParticlesAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Get the array of particles from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts = static_cast<AnaEventB*>(&event)->nParticles;

  // Count all particles with valid start positions
  box.nAllParticles = 0;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(!part) continue;

    // Check if particle has valid start position
    if (part->PositionStart[0] < -900 || part->PositionStart[1] < -900 || part->PositionStart[2] < -900) {
      continue; // Skip particles with invalid start positions
    }

    box.nAllParticles++;
  }

  return true;
}

//********************************************************************
bool HasEnoughParticlesCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have at least 2 particles with valid start positions
  return (box.nAllParticles >= 2);
}


//********************************************************************
bool FindNeutralCandidatesAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Clear existing neutral particle candidates
  box.neutralParticleCandidates.clear();
  box.nNeutralParticleCandidates = 0;
  box.nAnnihilationVerticesBeforeFiltering = 0;
  box.nAnnihilationVerticesAfterFiltering = 0;

  // Create annihilation vertices only in this branch.
  const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
  std::vector<AnaAnnihilationVertexPD*> annihilationVertices =
      pdAnnihilationUtils::CreateVertices(*static_cast<AnaEventB*>(&event), maxDaughterDistance,
                                          &box.nAnnihilationVerticesBeforeFiltering,
                                          &box.nAnnihilationVerticesAfterFiltering);

  // Build one wrapper candidate per annihilation vertex.
  box.neutralParticleCandidates =
      pdNeutralUtils::CreateAnnihilationOnlyNeutrals(*static_cast<AnaEventB*>(&event), annihilationVertices);
  box.nNeutralParticleCandidates = box.neutralParticleCandidates.size();

  return true;
}

//********************************************************************
bool HasNeutralCandidatesCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have at least one neutral particle candidate
  return (box.nNeutralParticleCandidates > 0);
}

//********************************************************************
bool K0StartEndDirCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0recostartenddir > 0.95
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    // Calculate dot product between start and end direction vectors
    Float_t k0recostartenddir = candidate->DirectionStart[0]*candidate->DirectionEnd[0] +
                                 candidate->DirectionStart[1]*candidate->DirectionEnd[1] +
                                 candidate->DirectionStart[2]*candidate->DirectionEnd[2];

    if (k0recostartenddir > 0.95) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

//********************************************************************
bool K0VtxOpeningCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0vtxrecoopening < 25
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    // Check if annihilation vertex exists
    if (!candidate->AnnihilationVertex) continue;

    Float_t k0vtxrecoopening = -999.0f;
    if (candidate->AnnihilationVertex->Particles.size() >= 2) {
      AnaParticlePD* d1 = candidate->AnnihilationVertex->Particles[0];
      AnaParticlePD* d2 = candidate->AnnihilationVertex->Particles[1];
      if (d1 && d2) {
        TVector3 dir1(d1->DirectionStart[0], d1->DirectionStart[1], d1->DirectionStart[2]);
        TVector3 dir2(d2->DirectionStart[0], d2->DirectionStart[1], d2->DirectionStart[2]);
        if (dir1.Mag() > 0 && dir2.Mag() > 0) {
          k0vtxrecoopening = static_cast<Float_t>(dir1.Angle(dir2));
        }
      }
    }

    // Check if opening angle is valid and passes cut
    if (k0vtxrecoopening > -900 && k0vtxrecoopening < 25.0) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

//********************************************************************
bool K0ParentTruePDGCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0partruepdg == 321
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    // Check if parent exists
    if (!candidate->Parent) continue;

    // Check if parent has true object
    if (!candidate->Parent->TrueObject) continue;

    // Get true parent PDG
    AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(candidate->Parent->TrueObject);
    if (!trueParent) continue;

    Int_t k0partruepdg = trueParent->PDG;

    if (k0partruepdg == 321) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

//********************************************************************
bool K0LengthCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0recolength > 5
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    Float_t k0recolength = candidate->Length;

    // Check if length is valid and passes cut
    if (k0recolength > -900 && k0recolength > 5.0) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

//********************************************************************
bool K0Dau1LengthCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0dau1recolength > 60
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    // Check if annihilation vertex exists and has at least one particle
    if (!candidate->AnnihilationVertex) continue;
    if (candidate->AnnihilationVertex->Particles.size() < 1) continue;

    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(candidate->AnnihilationVertex->Particles[0]);
    if (!daughter1) continue;

    Float_t k0dau1recolength = daughter1->Length;

    // Check if length is valid and passes cut
    if (k0dau1recolength > -900 && k0dau1recolength > 60.0) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

//********************************************************************
bool K0Dau2LengthCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have any candidates
  if (box.nNeutralParticleCandidates == 0) {
    return false;
  }

  // Check if at least one candidate passes the cut: k0dau2recolength < 180
  for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
    AnaNeutralParticlePD* candidate = box.neutralParticleCandidates[i];
    if (!candidate) continue;

    // Check if annihilation vertex exists and has at least two particles
    if (!candidate->AnnihilationVertex) continue;
    if (candidate->AnnihilationVertex->Particles.size() < 2) continue;

    AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(candidate->AnnihilationVertex->Particles[1]);
    if (!daughter2) continue;

    Float_t k0dau2recolength = daughter2->Length;

    // Check if length is valid and passes cut
    if (k0dau2recolength > -900 && k0dau2recolength < 180.0) {
      return true; // At least one candidate passes
    }
  }

  return false;
}

