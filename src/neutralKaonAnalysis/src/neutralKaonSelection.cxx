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

