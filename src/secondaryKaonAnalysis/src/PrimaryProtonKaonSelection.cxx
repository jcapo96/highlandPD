#include "PrimaryProtonKaonSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
PrimaryProtonKaonSelection::PrimaryProtonKaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void PrimaryProtonKaonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselection
  AddStep(StepBase::kAction,   "FIND PANDORA TRACK",    new FindBeamTrackAction(),        true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "BEAM TRACK EXISTS",     new BeamTrackExistsCut(),         true);// in pdBaseAnalysis/src/pandoraPreselection  
  //select events using beam instrumentation
  AddStep(StepBase::kCut,      "BEAM PDG FILTER",       new BeamPDGCut(2212),             true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "BEAM QUALITY",          new BeamQualityCut(),             true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kAction,   "GET PRIMARY PROTON",    new GetPrimaryProtonAction(),     true);
  AddStep(StepBase::kCut,      "CANDIDATE EXISTS",      new EventHasKaonCut(),            true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "DAUGHTER TRACK-LIKE",   new MuonIsTrackCut(),             true);
    AddStep(i, StepBase::kCut, "DAUGHTER CNN cut",      new MuonCNNCut(0.35,1.01),        true);
    AddStep(i, StepBase::kCut, "DAUGHTER #chi^{2}",     new MuonChi2Cut(0.5,6.00),        true);
    AddStep(i, StepBase::kCut, "DAUGHTER MOM",          new MuonRangeMomCut(0.221,0.245), true);
    AddStep(i, StepBase::kCut, "ANGLE",                 new KaonMuonAngleCut(-1,0.64),    true);
    AddStep(i, StepBase::kCut, "DISTANCE",              new KaonMuonDistanceCut(0,7.7),   true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    std::stringstream ssi;
    ssi << i;
    SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  }  

  // This number tells the selection to stop when a given accum level (accumulated cut level) is not passed. 
  // For example, 2 as argument would mean that selection is stopped for the current toy experiment when the 
  // third cut (numbering starts at 0) is not passed. This is a way of saving time.
  // -1 means no preselection
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool GetPrimaryProtonAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB); 
    
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part->isPandora && part->DaughtersIDs.size() == 1){
      box.Candidates.push_back(part);
      break; //only one pandora particle
    }
  }
  
  return true;  
}

//**************************************************
void PrimaryProtonKaonSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonCandidatesAndDaughters(event);
}

