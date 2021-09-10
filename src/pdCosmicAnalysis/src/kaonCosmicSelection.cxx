#include "kaonCosmicSelection.hxx"
#include "kaonCosmicAnalysis.hxx"
#include "kaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

/*
  This file contains:
  - the method DefineSteps in which actions and cuts are added to the selection in the right order. 
  - the implementation of each of the steps
  - the method InitializeEvent, in which the EventBox is created. This box is used to make a subsample 
    of event objects (Particles, hits, etc), which are frequently needed, such that we don't have to create 
    the subsample for every toy experiment dusing systematic propagation
 */


//********************************************************************
kaonCosmicSelection::kaonCosmicSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void kaonCosmicSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //secondary kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsAction()                  );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonFromCosmicCut(),  true);
  //split the selection in branches, one for each possible candidate
  AddSplit(kaonCosmicAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)kaonCosmicAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),             true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(),                true);
    AddStep(i, StepBase::kCut, "kaon daughter CNN cut",      new MuonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(),            true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon-muon angle cut",        new KaonMuonAngleCut(),           true);
    AddStep(i, StepBase::kCut, "kaon-muon distance cut",     new KaonMuonDistanceCut(),        true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)kaonCosmicAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
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
bool EventHasKaonFromCosmicCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //if the main track has at least a single daughter with a single daughter, pass
  if (box.Candidates.size()>0) return true;
  else return false;
}

//**************************************************
void kaonCosmicSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();
  
  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}
