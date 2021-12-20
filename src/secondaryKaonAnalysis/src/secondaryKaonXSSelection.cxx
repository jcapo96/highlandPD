#include "secondaryKaonXSSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

//********************************************************************
secondaryKaonXSSelection::secondaryKaonXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void secondaryKaonXSSelection::DefineSteps(){
//********************************************************************

  //copy steps from pandoraPreselection
  AddStep(StepBase::kAction,   "find Pandora track",         new FindPandoraTrackAction());// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "Pandora track exists",       new CandidateExistsCut()    );// in pdBaseAnalysis/src/pandoraPreselection  
  //select events using beam instrumentation
  AddStep(StepBase::kCut,      "beam pdg filter"     ,       new BeamFilterForXSCut()    );
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsForXSAction()   );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),           true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(10),            true);
    AddStep(i, StepBase::kCut, "kaon daughter CNN cut",      new MuonCNNCut(0.5),            true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(0.22,0.04), true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(0.5),            true);
    AddStep(i, StepBase::kCut, "kaon-muon angle cut",        new KaonMuonAngleCut(0.6),      true);
    AddStep(i, StepBase::kCut, "kaon-muon distance cut",     new KaonMuonDistanceCut(10),    true);
    AddStep(i, StepBase::kCut, "kaon chi2 cut",              new KaonChi2Cut(),              true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    std::stringstream ssi;
    ssi << i;
    SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  }

  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool BeamFilterForXSCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut* beamCut1 = new BeamPDGCut(211);
  BeamPDGCut* beamCut2 = new BeamPDGCut(13);

  if(beamCut1->Apply(event,boxB) || beamCut2->Apply(event,boxB))return true;
  else return false;
}

//**************************************************
bool GetKaonsForXSAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // Get the beam track
  AnaParticlePD* beamTrack = box.MainTrack;
    
  //look over daughters and look for candidates
  for(int i = 0; i < (int)beamTrack->Daughters.size(); i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(beamTrack->Daughters[i]);
    if(part->Daughters.size() == 1)box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool KaonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its chi2
  AnaParticlePD* kaon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]);
  if(kaon->Chi2Muon / kaon->Chi2ndf < 150)return true;
  else return false;
}

//**************************************************
void secondaryKaonXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}

