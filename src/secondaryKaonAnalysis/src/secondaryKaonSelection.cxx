#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryKaonSelection::secondaryKaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void secondaryKaonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselection
  AddStep(StepBase::kAction,   "FIND PANDORA TRACK",         new FindBeamTrackAction(),true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "BEAM TRACK EXISTS",           new BeamTrackExistsCut(), true);// in pdBaseAnalysis/src/pandoraPreselection  
  //select events using beam instrumentation
  //AddStep(StepBase::kCut,      "beam pdg filter",            new BeamHadronCut(),      true);
  AddStep(StepBase::kAction,   "GET VECTOR OF KAONS",      new GetKaonsAction(),     true);
  AddStep(StepBase::kCut,      "CANDIDATE EXISTS",             new EventHasKaonCut(),    true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "DAUGHTER TRACK-LIKE",   new MuonIsTrackCut(),             true);
    AddStep(i, StepBase::kCut, "DAUGHTER #chi^{2}",     new MuonChi2Cut(0.5,6.00),        true);
    AddStep(i, StepBase::kCut, "DAUGHTER MOM",          new MuonRangeMomCut(0.221,0.245), true);
    AddStep(i, StepBase::kCut, "ANGLE",                 new KaonMuonAngleCut(-1,0.64),    true);
    AddStep(i, StepBase::kCut, "DISTANCE",              new KaonMuonDistanceCut(0,7.7),   true);
    // AddStep(i, StepBase::kCut, "Proton chi2 cut",            new ProtonChi2Cut(0,84),          true);
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
bool BeamHadronCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut beamCut0(11);
  BeamPDGCut beamCut1(13);
  BeamPDGCut beamCut2(211);
  BeamPDGCut beamCut3(321);
  BeamPDGCut beamCut4(2212);

  if(beamCut0.Apply(event,boxB) ||
     beamCut1.Apply(event,boxB) || beamCut2.Apply(event,boxB) || 
     beamCut4.Apply(event,boxB) || beamCut3.Apply(event,boxB))
    return true;
  else return false;
}

//**************************************************
bool GetKaonsAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
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
    if(part->isPandora)continue; //skip beam particle
    if(part->DaughtersIDs.size() == 1 && part->ParentID != -1)box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool EventHasKaonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //if the main track has at least a single daughter with a single daughter, pass
  if (box.Candidates.size()>0) return true;
  else return false;
}

//**************************************************
bool MuonIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();
  
  //make sure there is a candidate for this branch
  if(branchesIDs[0] > box.Candidates.size()-1)return false;

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(!dau)return false;
  if(dau->Type == 2){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool MuonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
 
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter chi2
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->Chi2Muon / dau->Chi2ndf > _lower_cut && dau->Chi2Muon / dau->Chi2ndf < _upper_cut && dau->Chi2Muon > 0){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool MuonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter CNN score
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->CNNscore[0] > _lower_cut && dau->CNNscore[0] < _upper_cut){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool MuonRangeMomCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter range mom
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(mom > _lower_cut && mom < _upper_cut){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool KaonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on its CNN
  if(box.Candidates[branchesIDs[0]]->CNNscore[0] > _lower_cut){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool KaonMuonAngleCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on the angle between it and its daughter muon
  AnaParticlePD* kaon = box.Candidates[branchesIDs[0]];
  AnaParticlePD* muon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + kaon->DirectionEnd[i] * muon->DirectionStart[i];
  if(cos > _lower_cut && cos < _upper_cut){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool KaonMuonDistanceCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on the distance between it and its daughter muon
  AnaParticlePD* kaon = box.Candidates[branchesIDs[0]];
  AnaParticlePD* muon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(kaon->PositionEnd[i] - muon->PositionStart[i],2);
  dis = sqrt(dis);
  if(dis > _lower_cut && dis < _upper_cut){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false; 
}

//**************************************************
bool ProtonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //make sure there is a candidate for this branch
  if(branchesIDs[0] > box.Candidates.size()-1)return false;

  //if there is kaon, cut in its chi2
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]);
  double chi2 = part->Chi2Proton;
  double ndf  = part->Chi2ndf;
  if(chi2/ndf > _lower_cut && chi2/ndf < _upper_cut && chi2 > 0 && chi2>0){
    box.UpdateBestCandidateIndex(Index(),branchesIDs[0]);
    return true;
  }
  else return false;
}

//**************************************************
void secondaryKaonSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonCandidatesAndDaughters(event);
}

