#include "secondaryKaonSelectionNoBranches.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryKaonSelectionNoBranches::secondaryKaonSelectionNoBranches(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void secondaryKaonSelectionNoBranches::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselection
  AddStep(StepBase::kAction,   "FIND PANDORA TRACK",         new FindBeamTrackAction(),true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "BEAM TRACK EXISTS",          new BeamTrackExistsCut(), true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kAction,   "GET VECTOR OF KAONS",        new GetKaonsAction(),     true);
  AddStep(StepBase::kCut,      "CANDIDATE EXISTS",           new EventHasKaonCut(),    true);
  AddStep(StepBase::kCut,      "EVENT HAS SUPER CANDIDATE",  new SuperCandidateCut(),  true);

  SetBranchAlias(0,"NO BRANCH :)");

  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool SuperCandidateCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  bool SuperCandidate = false;

  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //loop over candidates
  for(int i = 0; i < (int)box.Candidates.size(); i++){
    if(secondaryKaonSelectionNoBranches::IsSuperCandidate(box.Candidates[i])){
      SuperCandidate = true;
      box.Candidates.erase(box.Candidates.begin()+1,box.Candidates.end());
    }
    else{
      box.Candidates.erase(box.Candidates.begin()+i);
      i--;
    }
  }

  box.Candidates.shrink_to_fit();
  box.UpdateBestCandidateIndex(Index(),0);

  return SuperCandidate;
}

//**************************************************
bool secondaryKaonSelectionNoBranches::IsSuperCandidate(AnaParticlePD* part){
//**************************************************

  return (MuonIsTrack(part) && MuonChi2(part) && MuonMom(part) &&
	  KaonMuonAngle(part) && KaonMuonDistance(part));
}

//**************************************************
bool secondaryKaonSelectionNoBranches::MuonIsTrack(AnaParticlePD* part){
//**************************************************

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return false;
  if(dau->Type == 2)
    return true;
  else return false; 
}

//**************************************************
bool secondaryKaonSelectionNoBranches::MuonChi2(AnaParticlePD* part){
//**************************************************

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return false;
  if(dau->Chi2Muon / dau->Chi2ndf > 0.5 && dau->Chi2Muon / dau->Chi2ndf < 6.0 && dau->Chi2Muon > 0)
    return true;
  else return false; 
}

//**************************************************
bool secondaryKaonSelectionNoBranches::MuonMom(AnaParticlePD* part){
//**************************************************

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return false;
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(mom > 0.221 && mom < 0.245)
    return true;
  else return false; 
}

//**************************************************
bool secondaryKaonSelectionNoBranches::KaonMuonAngle(AnaParticlePD* part){
//**************************************************
  
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau || !part)return false;
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + part->DirectionEnd[i] * dau->DirectionStart[i];
  if(cos > -1 && cos < 0.64)
    return true;  
  else return false; 
}

//**************************************************
bool secondaryKaonSelectionNoBranches::KaonMuonDistance(AnaParticlePD* part){
//**************************************************
  
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau || !part)return false;
  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(part->PositionEnd[i] - dau->PositionStart[i],2);
  dis = sqrt(dis);
  if(dis > 0 && dis < 7.7)
    return true;  
  else return false; 
}

//**************************************************
void secondaryKaonSelectionNoBranches::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonCandidatesAndDaughters(event);
}

