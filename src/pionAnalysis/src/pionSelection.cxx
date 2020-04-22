#include "pionSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pionAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

//********************************************************************
pionSelection::pionSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void pionSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselecion
  AddStep(StepBase::kAction, "find Pandora track",         new FindPandoraTrackAction());  
  AddStep(StepBase::kCut,    "candidate exists",           new CandidateExistsCut());
  AddStep(StepBase::kCut,    "beam is pion",               new BeamPionCut());
  AddStep(StepBase::kCut,    "pandora reco worked",        new CandidateIsBeamCut());
  AddStep(StepBase::kCut,    "candidate is track",         new CandidateIsTrackCut());
  AddStep(StepBase::kCut,    "pion track end",             new PionEndsAPA3Cut());
  AddStep(StepBase::kCut,    "seltrk chi2 cut",            new PionPassChi2Cut());
  AddStep(StepBase::kAction, "compute daughter distance",  new ComputeDaughterDistanceAction());  
  AddStep(StepBase::kCut,    "no pion daughter",           new NoPionDaughterCut());

  //Add a split to the trunk with 2 branches.
  AddSplit(2);

  //First branch is for pi0 showers
  AddStep(0, StepBase::kCut, "pi0 showers",        new PionCexCut());
  //Second branch is for pi absortion
  AddStep(1, StepBase::kCut, "pi absorption",      new PionAbsorptionCut());

  // Set the branch aliases to the three branches
  SetBranchAlias(0,"pi0 shower",     0);
  SetBranchAlias(1,"pi absorption",  1);

  // No preselection for the moment
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool BeamPionCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)boxB;
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  if (trueBeamPart){
    // for MC
    if (trueBeamPart->PDG==211 || trueBeamPart->PDG==-13) return true;
  }
  else{
    // for real DATA
    for(int i = 0; i < (int)beam->PDGs.size(); i++){
      if (beam->PDGs[i] == 211) return true;
    }
  }
  return false;
}

//**************************************************
bool CandidateIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the maintrack
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.MainTrack);
  
  //if the seltrk is a track, accept it
  if (part->Type == AnaParticlePD::kTrack) return true;
  else return false;
}

//**************************************************
bool PionEndsAPA3Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cutAPA3_Z = 226.;

  // Cast the ToyBox to the appropiate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  if(box.MainTrack->PositionEnd[2] < cutAPA3_Z)return true;
  else return false;  
}

//**************************************************
bool PionPassChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cut_chi2 = 140;

  //Cast the ToyBox to the appropiate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);
  if (!box.MainTrack) return false;
  
  if(box.MainTrack->Chi2Proton/box.MainTrack->Chi2ndf > cut_chi2) return true;
  else return false;
}

//**************************************************
bool ComputeDaughterDistanceAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 
    
  //compute distance between daughters and vertex
  pdAnaUtils::ComputeDistanceToVertex(box.MainTrack, box.DaughterDistanceToVertex);

  return true;  
}

//**************************************************
bool NoPionDaughterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cut_daughter_track_distance = 10;
  Float_t cut_CNNTrackScore = 0.35;
  Float_t cut_chi2 = 50;
  
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  // Look for pion daughters (track hypothesis)
  bool noPion = true;
  for(UInt_t i = 0; i < box.MainTrack->Daughters.size()/2; i++){
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);
    if(daughter->UniqueID != -999 &&
       daughter->CNNscore[0]                  > cut_CNNTrackScore &&
       daughter->Chi2Proton/daughter->Chi2ndf > cut_chi2 &&
       box.DaughterDistanceToVertex[i]        < cut_daughter_track_distance){
      noPion = false;
      break;
    }
  }
  return noPion;  
}

//**************************************************
bool PionCexCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  Float_t cut_daughter_shower_distance_low  = 2.;
  Float_t cut_daughter_shower_distance_high = 100.;
  int cut_nHits_shower_low  = 12;
  int cut_nHits_shower_high = 1000;
  Float_t cut_CNNTrackScore = 0.35;
  
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  // loop over daughters assuming shower hypothesis
  bool cex = false;
  for(UInt_t i = box.MainTrack->Daughters.size()/2; i < box.MainTrack->Daughters.size(); i++){
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);
    if(daughter->UniqueID != -999 &&
       daughter->CNNscore[0]           < cut_CNNTrackScore &&
       box.DaughterDistanceToVertex[i] > cut_daughter_shower_distance_low &&
       box.DaughterDistanceToVertex[i] < cut_daughter_shower_distance_high &&
       daughter->NHits                 > cut_nHits_shower_low &&
       daughter->NHits                 < cut_nHits_shower_high){
      cex = true;
      break;
    }
  }
  return cex;
}

//**************************************************
bool PionAbsorptionCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  PionCexCut cut;
  return !cut.Apply(event, boxB);
}

//**************************************************
void pionSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
}
