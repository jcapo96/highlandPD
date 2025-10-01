#include "neutralKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"
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

  // First create reconstructed vertices
  const double maxVertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");
  const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
  std::vector<AnaVertexPD*> vertices = pdAnaUtils::CreateReconstructedVertices(*static_cast<AnaEventB*>(&event), maxVertexRadius, maxDaughterDistance);

  // Then create neutral particle candidates from the vertices
  const double vertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");
  const double impactParameter = ND::params().GetParameterD("neutralKaonAnalysis.ImpactParameter");
  box.neutralParticleCandidates = pdAnaUtils::CreateAnaNeutralParticles(*static_cast<AnaEventB*>(&event), vertices, vertexRadius, impactParameter);
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

