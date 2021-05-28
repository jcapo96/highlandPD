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
kaonSelection::kaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void kaonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselecion
  AddStep(StepBase::kAction,   "find Pandora track",         new FindPandoraTrackAction());// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "candidate exists",           new CandidateExistsCut()    );// in pdBaseAnalysis/src/pandoraPreselection  
  //select events using beam instrumentation
  AddStep(StepBase::kCut,      "beam pdg filter",            new BeamFilterCut()         );
  //secondary kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsAction()        );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  //AddStep(StepBase::kAction,   "get a vector of daughters",  new GetDaughtersAction()    );
  //split the selection in branches, one for each possible candidate
  AddSplit(kaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)kaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    //AddStep(i, StepBase::kCut, "daughter is kaon candidate", new DaughterIsKaonCandidateCut(), true);
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),             true);
    //AddStep(i, StepBase::kCut, "kaon daughter is muon-like", new MuonFromDecayCut(),           true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(),                true);
    AddStep(i, StepBase::kCut, "kaon daughter CNN cut",      new MuonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(),            true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon-muon angle cut",        new KaonMuonAngleCut(),           true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)kaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
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
bool BeamFilterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)boxB;
  
  //get the beam
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  //for MC
  if(event.GetIsMC()){
    if(beam->BeamParticle){
      if(beam->BeamParticle->TrueObject){
	//cast the beam particle
	AnaTrueParticlePD* beamTruePart = static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);
        if(beamTruePart->PDG==211 || beamTruePart->PDG==321 || beamTruePart->PDG==2212)return true;
      }
    }
  }
  
  return false;
}


//**************************************************
bool GetKaonsAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 
    
  // Loop over daughters and granddaughters add them to the candidates vector
  // if they have a single daughter
  for(int i = 0; i < (int)box.MainTrack->Daughters.size(); i++){
    //cast the daughter in the appropiate class
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);
    if(dau->Daughters.size() == 1)box.Candidates.push_back(dau);
    for(int j = 0; j < (int)dau->Daughters.size(); j++){
      //cast the daughter in the appropiate class
      AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[j]);
      if(gdau->Daughters.size() == 1)box.Candidates.push_back(gdau);
    }
  }
  
  return true;  
}

//**************************************************
bool GetDaughtersAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 
    
  // Loop over daughters and add them to the candidates vector
  for(int i = 0; i < std::min((int)box.MainTrack->Daughters.size(),(int)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS); i++){
    //cast the daughter in the appropiate class
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);
    box.Candidates.push_back(dau);
  }
  
  return true;  
}


//**************************************************
bool EventHasKaonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  // Main track must exist
  if (!box.MainTrack) return false;
  
  //if the main track has at least a single daughter with a single daughter, pass
  if (box.Candidates.size()>0) return true;
  else return false;
}

//**************************************************
bool DaughterIsKaonCandidateCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  // Main track must exist
  if (!box.MainTrack) return false;
  
  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //check if there is no daughter for this branch
  if(branchesIDs[0] >= box.Candidates.size())return false;

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]);
  if(dau->Daughters.size() == 1)return true;
  else return false; 
}

//**************************************************
bool MuonIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();
  
  //make sure there is a candidate for this branch
  if(branchesIDs[0] > box.Candidates.size()-1)return false;

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->Type == 2)return true;
  else return false; 
}

//**************************************************
bool MuonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter chi2
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->Chi2Muon / dau->Chi2ndf < 10)return true;
  else return false; 
}

//**************************************************
bool MuonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter CNN score
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->CNNscore[0] > 0.6)return true;
  else return false; 
}

//**************************************************
bool MuonRangeMomCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter range mom
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(abs(mom-0.22) < 0.02)return true;
  else return false; 
}

//**************************************************
bool MuonFromDecayCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, check that its daughter is a muon from kaon decay
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(dau->Chi2Muon / dau->Chi2ndf < 10 &&
     dau->CNNscore[0] > 0.6 && 
     abs(mom-0.22) < 0.02)return true;
  else return false; 
}

//**************************************************
bool KaonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on its CNN
  if(box.Candidates[branchesIDs[0]]->CNNscore[0] > 0.8)return true;
  else return false; 
}

//**************************************************
bool KaonMuonAngleCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on the angle between it and its daughter muon
  AnaParticlePD* kaon = box.Candidates[branchesIDs[0]];
  AnaParticlePD* muon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + kaon->DirectionEnd[i] * muon->DirectionStart[i];
  if(cos < 0.1 )return true;
  else return false; 
}


//**************************************************
void kaonSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}

