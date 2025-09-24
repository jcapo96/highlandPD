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
  AddStep(StepBase::kCut, "beam pdg filter", new BeamPDGCut(321), true);

  // Preliminary K0 Selection Steps - Look for vertex candidates
  AddStep(StepBase::kAction, "find all particles"   , new FindAllParticlesAction(), true);
  AddStep(StepBase::kCut   , "enough particles" , new HasEnoughParticlesCut(), true);
  AddStep(StepBase::kAction, "find vertex candidates", new FindVertexCandidatesAction(), true);
  AddStep(StepBase::kCut   , "has vertex candidates" , new HasVertexCandidatesCut(), true);

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
bool GetNeutralKaonsAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
  //********************************************************************
  (void)event;

  // Cast the box
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts = static_cast<AnaEventB*>(&event)->nParticles;

  if (!parts || nParts <= 0) {
      std::cerr << "Error: No particles found in event!" << std::endl;
      return false;
  }

  // Step 1: Find all Kaons and store in Candidates[0]
  box.ResetBase();
  for (int i = 0; i < nParts; i++) {
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if (!part || part->DaughtersIDs.empty()) continue;  // Validate pointer
  }
  return true;
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
bool FindVertexCandidatesAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Clear existing vertex candidates
  box.trueVertexCandidates.clear();
  box.reconVertexCandidates.clear();
  box.nTrueVertexCandidates = 0;
  box.nReconVertexCandidates = 0;

  // Use the common utility function to create reconstructed vertices
  const double maxVertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");
  const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
  box.reconVertexCandidates = pdAnaUtils::CreateReconstructedVertices(*static_cast<AnaEventB*>(&event), maxVertexRadius, maxDaughterDistance);
  box.nReconVertexCandidates = box.reconVertexCandidates.size();

  // Create true vertex candidates from the reconstructed vertices
  for (const auto& reconVertex : box.reconVertexCandidates) {
    if (!reconVertex || reconVertex->NParticles != 3) continue; // Skip if not 1 parent + 2 daughters

    // Create true vertex candidate if all particles have true information
    std::vector<AnaTrueParticlePD*> trueVertexParticles;
    bool allHaveTrueInfo = true;

    for (int i = 0; i < reconVertex->NParticles; i++) {
      AnaParticlePD* particle = static_cast<AnaParticlePD*>(reconVertex->Particles[i]);
      if (!particle) {
        allHaveTrueInfo = false;
        break;
      }

      AnaTrueParticleB* truePart = particle->GetTrueParticle();
      if (truePart) {
        trueVertexParticles.push_back(static_cast<AnaTrueParticlePD*>(truePart));
      } else {
        allHaveTrueInfo = false;
        break;
      }
    }

    if (allHaveTrueInfo && trueVertexParticles.size() == 3) { // exactly 1 parent + 2 daughters
      AnaTrueVertexPD* trueVertex = new AnaTrueVertexPD();

      // Get the true parent particle (first particle in the vertex)
      AnaTrueParticleB* trueParent = static_cast<AnaParticlePD*>(reconVertex->Particles[0])->GetTrueParticle();
      if (trueParent) {
        // Set true vertex position to the true parent's position
        trueVertex->Position[0] = trueParent->Position[0];
        trueVertex->Position[1] = trueParent->Position[1];
        trueVertex->Position[2] = trueParent->Position[2];
        trueVertex->Position[3] = 0.0; // time

        // Set true vertex properties - exactly 1 parent + 2 daughters
        trueVertex->NTrueParticles = 3; // parent + 2 daughters
        trueVertex->TrueParticles = trueVertexParticles;
        trueVertex->Parent = static_cast<AnaTrueParticlePD*>(trueParent);
        trueVertex->Generation = 1; // Secondary vertex
        trueVertex->ReactionType = 0; // Unknown for now

        // Add to true candidates
        box.trueVertexCandidates.push_back(trueVertex);
        box.nTrueVertexCandidates++;
      }
    }
  }


  return true;
}

//********************************************************************
bool HasVertexCandidatesCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have at least one vertex candidate (either reconstructed or true)
  return (box.nReconVertexCandidates > 0 || box.nTrueVertexCandidates > 0);
}

//********************************************************************
VertexDaughterCountCut::VertexDaughterCountCut(int min_daughters) :
  _min_daughters(min_daughters) {
//********************************************************************
}

//********************************************************************
bool VertexDaughterCountCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check reconstructed vertex candidates
  for (const auto& reconVertex : box.reconVertexCandidates) {
    if (reconVertex) {
      // Count daughters (particles that are not the parent)
      int daughterCount = 0;
      for (int i = 0; i < reconVertex->NParticles; i++) {
        if (reconVertex->Particles[i] && reconVertex->Particles[i] != reconVertex->Parent) {
          daughterCount++;
        }
      }

      if (daughterCount >= _min_daughters) {
        return true; // Found at least one valid vertex
      }
    }
  }

  // Check true vertex candidates
  for (const auto& trueVertex : box.trueVertexCandidates) {
    if (trueVertex) {
      // For true vertices, we can't easily count daughters separately
      // This is a limitation of the current data structure
      // For now, we'll skip this check for true vertices
      continue;
    }
  }

  return false; // No valid vertices found
}

//********************************************************************
VtxDaughtersAreDaughtersOfVtxParentCut::VtxDaughtersAreDaughtersOfVtxParentCut() {
//********************************************************************
}

//********************************************************************
bool VtxDaughtersAreDaughtersOfVtxParentCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Get the current branch
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  // Make sure there is a candidate for this branch
  if(branchesIDs[0] > box.reconVertexCandidates.size()-1) return false;

  // Get the vertex candidate for this branch
  AnaVertexPD* reconVertex = box.reconVertexCandidates[branchesIDs[0]];
  if (!reconVertex || reconVertex->NParticles < 3) return false;

  // Check if all daughters are daughters of the parent
  AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(reconVertex->Particles[1]->TrueObject);
  AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(reconVertex->Particles[2]->TrueObject);
  AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(reconVertex->Particles[0]->TrueObject);

  if (trueDaughter1 && trueDaughter2 && trueParent) {
    if (trueDaughter1->ParentID == trueParent->ID && trueDaughter2->ParentID == trueParent->ID) {
      box.UpdateBestCandidateIndex(Index(), branchesIDs[0]);
      return true; // Found a valid vertex for this branch
    }
  }

  return false; // No valid vertex found for this branch
}



