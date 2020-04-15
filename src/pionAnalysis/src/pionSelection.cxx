#include "pionSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pionAnalysisUtils.hxx"

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
  AddStep(StepBase::kAction, "find main track",            new FindBeamPionTrackAction());  
  AddStep(StepBase::kCut,    "beam pion",                  new BeamPionCut());
  AddStep(StepBase::kCut,    "beam is track",              new BeamIsTrackCut());
  AddStep(StepBase::kCut,    "beam geometric",             new BeamPionGeometricCut());
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
bool FindBeamPionTrackAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 

  // This action fills box.MainTrack 
  box.MainTrack = NULL;
  
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  
  // loop ever particles in the event
  for (Int_t i=0;i<nParts; ++i){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //save the ones that are correctly selected by Pandora
    bool beamtrack=false;
    if (useIsBeamLike)      
      beamtrack = pdAnaUtils::isBeamLike(part,beam);
    else
      beamtrack = part->isBeamPart;

    if(beamtrack){
      box.MainTrack = part;
      break;
    }
  }  
  return true;
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
      if (beam->PDGs[i] == 211)return true;
    }
  }
  return false;
}

//**************************************************
bool BeamIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  //if the main track exists cast it
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.MainTrack);
  
  //if the seltrk is a track, accept it
  if (part->Type == AnaParticlePD::kTrack) return true;
  else return false;
}

//**************************************************
bool BeamPionGeometricCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  //From Owen Goodwins studies
  Float_t mccuts[7]  ={-3.,  7., -8.,  7., 27.5, 32.5, 0.93};
  Float_t datacuts[7]={ 0., 10., -5., 10., 30.,  35.0, 0.93};
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  //if the beam part exists cast it
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  if (!beam->BeamParticle) return false;
  AnaParticle* beampart = static_cast<AnaParticle*>(beam->BeamParticle);

  Float_t beampos[3],beamdir[3], dist[3], dcos=0, cuts[7];      

  // different way of obtaining the beam position and angle for DATA and MC
  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beampart->TrueObject);

  // For MC
  if (trueBeamPart){
    for (int i=0;i<3;i++){
      beampos[i] = trueBeamPart->Position[i]-trueBeamPart->Position[2]*(trueBeamPart->Direction[i]/trueBeamPart->Direction[2]);
      beamdir[i] = trueBeamPart->Direction[i];
    }
    for (int i=0;i<7;i++) cuts[i] = mccuts[i];
  }
  else{
    // For Data
    if(beam->nMomenta != 1 || beam->nTracks != 1)return false;
    for (int i=0;i<3;i++){
      beampos[i] = beampart->PositionEnd[i];
      beamdir[i] = beampart->DirectionEnd[i];
    }
    for (int i=0;i<7;i++) cuts[i] = datacuts[i];
  }

  // compute the difference in position and cos(angle)
  for (int i=0;i<3;i++){
    dist[i] = box.MainTrack->PositionStart[i] - beampos[i];
    dcos   += box.MainTrack->DirectionStart[i]*beamdir[i];
  }
  
  if(dist[0] < cuts[0] || dist[0] > cuts[1] || dist[1] < cuts[2] || dist[1] > cuts[3] || dist[2] < cuts[4] || dist[2] > cuts[5] || dcos < cuts[6]) return false;
  else return true;
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
    AnaParticle* daughter = static_cast<AnaParticle*>(box.MainTrack->Daughters[i]);
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
    AnaParticle* daughter = static_cast<AnaParticle*>(box.MainTrack->Daughters[i]);
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
