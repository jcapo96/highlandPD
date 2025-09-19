#include "neutralKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
neutralKaonSelection::neutralKaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void neutralKaonSelection::DefineSteps(){
//********************************************************************

  //basic cuts → see pdBaseAnalysis/src/pdBaseSelection
  AddStep(StepBase::kAction, "find Pandora track"    , new FindBeamTrackAction(), true);
  AddStep(StepBase::kCut   , "beam track in TPC"     , new BeamTrackExistsCut() , true);

  // Preliminary K0 Selection Steps
  AddStep(StepBase::kCut, "find beam daughters"   , new FindBeamDaughtersAction(), true);
  AddStep(StepBase::kCut   , ">2 daughters" , new HasEnoughBeamDaughtersCut(), true);
  AddStep(StepBase::kAction, "check K0 in truth"     , new CheckK0InTruthAction(), true);
  AddStep(StepBase::kCut   , ">3cm from end", new BeamDaughtersDistanceCut(3.0, 50.0), true);
  AddStep(StepBase::kCut   , "close pairs"    , new BeamDaughtersPairDistanceCut(0.0, 1.5), true);

  // Original neutral kaon steps (commented out for now)
  // AddStep(StepBase::kAction,   "Get Triplet of (parent, daughter1, daughter2)",      new GetNeutralKaonsAction(),     true);
  // AddStep(StepBase::kCut   , "beam pdg filter KPos"  , new BeamPDGCut(321)      , true);
  // AddStep(StepBase::kCut   , "BEAM quality cut"      , new BeamQualityCut()     , true);

  SetBranchAlias(0,"trunk");
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
  // box.neutralKaonCandidates.clear();
  // box.Candidates.shrink_to_fit();
  box.ResetBase();
  std::cout << "Size of Candidates before clear: " << box.neutralKaonCandidates.size() << std::endl;
  for (int i = 0; i < nParts; i++) {
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if (!part || part->DaughtersIDs.empty()) continue;  // Validate pointer

    std::cout << "Found a Kaon" << std::endl;
    std::cout << "Number of particles: " << nParts << std::endl;
    box.neutralKaonCandidates.push_back(part);
  }
  return true;
}

//********************************************************************
bool FindBeamDaughtersAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Get the array of particles from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts = static_cast<AnaEventB*>(&event)->nParticles;

  // Get beam particle
  AnaEventB& eventB = *static_cast<AnaEventB*>(&event);
  if (!eventB.Beam) return false;

  AnaBeam* beam = static_cast<AnaBeam*>(eventB.Beam);
  if (!beam || !beam->BeamParticle) return false;

  AnaParticleB* beamParticle = beam->BeamParticle;
  AnaTrueParticleB* beamTrue = beamParticle->GetTrueParticle();
  if (!beamTrue) return false;

  // Look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(!part) continue;

    // Get the true particle information
    AnaTrueParticleB* truePart = part->GetTrueParticle();
    if(!truePart) continue;

    // Check if this reconstructed particle is a daughter of the beam particle
    if(truePart->ParentID == beamTrue->ID){
      box.nBeamDaughters++;
    }
  }

  return true;
}

//********************************************************************
bool CheckK0InTruthAction::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Get beam particle
  AnaEventB& eventB = *static_cast<AnaEventB*>(&event);
  if (!eventB.Beam) return false;

  AnaBeam* beam = static_cast<AnaBeam*>(eventB.Beam);
  if (!beam || !beam->BeamParticle) return false;

  AnaParticleB* beamParticle = beam->BeamParticle;
  AnaTrueParticleB* beamTrue = beamParticle->GetTrueParticle();
  if (!beamTrue) return false;

  // Check for K0 in truth information that is a daughter of beam particle
  box.hasK0InTruth = false;
  AnaTrueParticleB** trueParticles = static_cast<AnaEventB*>(&event)->TrueParticles;
  int nTrueParticles = static_cast<AnaEventB*>(&event)->nTrueParticles;

  for (int i = 0; i < nTrueParticles; i++) {
    AnaTrueParticleB* truePart = trueParticles[i];
    if (truePart && beamTrue && beamTrue->ID >= 0 &&
        truePart->ParentID == beamTrue->ID && truePart->PDG == 310) {
      box.hasK0InTruth = true;
      break;
    }
  }

  return true;
}

//********************************************************************
bool HasEnoughBeamDaughtersCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;

  // Cast the box to the appropriate type
  ToyBoxNeutralKaon& box = *static_cast<ToyBoxNeutralKaon*>(&boxB);

  // Check if we have at least 2 beam daughters
  return (box.nBeamDaughters >= 2);
}

//********************************************************************
bool BeamDaughtersDistanceCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;
  (void)boxB; // Unused but required for framework

  // Get beam particle
  AnaEventB& eventB = *static_cast<AnaEventB*>(&event);
  if (!eventB.Beam) return false;

  AnaBeam* beam = static_cast<AnaBeam*>(eventB.Beam);
  if (!beam || !beam->BeamParticle) return false;

  AnaParticleB* beamParticle = beam->BeamParticle;
  AnaTrueParticleB* beamTrue = beamParticle->GetTrueParticle();
  if (!beamTrue) return false;

  // Get beam particle end position
  float beamEndX = beamTrue->PositionEnd[0];
  float beamEndY = beamTrue->PositionEnd[1];
  float beamEndZ = beamTrue->PositionEnd[2];

  // Check if beam particle end position is valid
  if (beamEndX < -900 || beamEndY < -900 || beamEndZ < -900) {
    return false; // Skip events with invalid beam particle end position
  }

  // Get the array of particles from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts = static_cast<AnaEventB*>(&event)->nParticles;

  // Find daughters of the beam particle and check their distance from beam end
  int validDaughters = 0;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(!part) continue;

    // Get the true particle information
    AnaTrueParticleB* truePart = part->GetTrueParticle();
    if(!truePart) continue;

    // Check if this reconstructed particle is a daughter of the beam particle
    if(truePart->ParentID == beamTrue->ID){
      // Check if daughter has valid start position
      if (part->PositionStart[0] < -900 || part->PositionStart[1] < -900 || part->PositionStart[2] < -900) {
        continue; // Skip daughters with invalid start positions
      }

      // Calculate distance from beam particle end to daughter start position
      double dist = sqrt((part->PositionStart[0] - beamEndX) * (part->PositionStart[0] - beamEndX) +
                         (part->PositionStart[1] - beamEndY) * (part->PositionStart[1] - beamEndY) +
                         (part->PositionStart[2] - beamEndZ) * (part->PositionStart[2] - beamEndZ));

      // Check if daughter is within the distance range from beam end
      if (dist >= _lower_cut && dist <= _upper_cut) {
        validDaughters++;
      }
    }
  }

  // Need at least 2 valid daughters within the distance range
  return (validDaughters >= 2);
}

//********************************************************************
bool BeamDaughtersPairDistanceCut::Apply(AnaEventC& event, ToyBoxB& boxB) const {
//********************************************************************

  (void)event;
  (void)boxB; // Unused but required for framework

  // Get beam particle
  AnaEventB& eventB = *static_cast<AnaEventB*>(&event);
  if (!eventB.Beam) return false;

  AnaBeam* beam = static_cast<AnaBeam*>(eventB.Beam);
  if (!beam || !beam->BeamParticle) return false;

  AnaParticleB* beamParticle = beam->BeamParticle;
  AnaTrueParticleB* beamTrue = beamParticle->GetTrueParticle();
  if (!beamTrue) return false;

  // Get the array of particles from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts = static_cast<AnaEventB*>(&event)->nParticles;

  // Find all daughters of the beam particle
  std::vector<AnaParticlePD*> beamDaughters;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(!part) continue;

    // Get the true particle information
    AnaTrueParticleB* truePart = part->GetTrueParticle();
    if(!truePart) continue;

    // Check if this reconstructed particle is a daughter of the beam particle
    if(truePart->ParentID == beamTrue->ID){
      // Check if daughter has valid start position
      if (part->PositionStart[0] < -900 || part->PositionStart[1] < -900 || part->PositionStart[2] < -900) {
        continue; // Skip daughters with invalid start positions
      }
      beamDaughters.push_back(part);
    }
  }

  // Need at least 2 daughters to check pair distance
  if (beamDaughters.size() < 2) return false;

  // Check all possible pairs of daughters for distance between their start positions
  for (UInt_t d1 = 0; d1 < beamDaughters.size(); d1++) {
    for (UInt_t d2 = d1 + 1; d2 < beamDaughters.size(); d2++) {
      AnaParticlePD* dau1 = beamDaughters[d1];
      AnaParticlePD* dau2 = beamDaughters[d2];

      if (!dau1 || !dau2) continue;

      // Calculate distance between the two daughters' start positions
      double distBetweenDaughters = sqrt((dau1->PositionStart[0] - dau2->PositionStart[0]) * (dau1->PositionStart[0] - dau2->PositionStart[0]) +
                                         (dau1->PositionStart[1] - dau2->PositionStart[1]) * (dau1->PositionStart[1] - dau2->PositionStart[1]) +
                                         (dau1->PositionStart[2] - dau2->PositionStart[2]) * (dau1->PositionStart[2] - dau2->PositionStart[2]));

      // Check if daughters start positions are close (within 3.0 cm as per PreliminaryK0Selection.C)
      if (distBetweenDaughters >= _lower_cut && distBetweenDaughters <= _upper_cut) {
        return true; // Found at least one valid pair
      }
    }
  }

  return false; // No valid pairs found
}

