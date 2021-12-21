#include "stoppingMuonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
stoppingMuonSelection::stoppingMuonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void stoppingMuonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")
  AddStep(StepBase::kAction, "find Pandora track",  new FindBeamTrackAction());
  AddStep(StepBase::kCut,    "beam track in TPC",   new BeamTrackExistsCut());
  AddStep(StepBase::kCut,    "beam muon",           new BeamMuonCut());
  //AddStep(StepBase::kCut,    "pandora reco worked", new CandidateIsBeamCut());
  AddStep(StepBase::kCut,    "muon CSDA range",     new MuonCSDARangeCut());
  
  SetBranchAlias(0,"trunk");
}

//**************************************************
bool BeamMuonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  BeamPDGCut* beamCut1 = new BeamPDGCut(211);
  BeamPDGCut* beamCut2 = new BeamPDGCut(13);

  if(beamCut1->Apply(event,boxB) || beamCut2->Apply(event,boxB))return true;
  else return false;
}

//**************************************************
bool BeamMuonAngleCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  //Cast the particle and define vector
  AnaParticle* part = static_cast<AnaParticle*>(box.MainTrack);
  double beamdir[3] = {-0.18,-0.20,0.96};//beam direction
  double partdir[3] = {part->DirectionStart[0],part->DirectionStart[1],part->DirectionStart[2]};

  //Cut
  if((beamdir[0]*partdir[0]+beamdir[1]*partdir[1]+beamdir[2]*partdir[2])>0.9)return true;
  else return false;  
}

//**************************************************
bool MuonCSDARangeCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  //check if it exists a beam particle
  AnaBeam* beam         = static_cast<AnaBeam*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;
  if(beamPart){
    Float_t mom = beamPart->Momentum;
    if (beamPart->TrueObject) mom =  static_cast<AnaTrueParticle*>(beamPart->TrueObject)->Momentum;
    Float_t length = static_cast<AnaParticlePD*>(box.MainTrack)->Length;
    Float_t csdarange = pdAnaUtils::ComputeCSDARange(mom*1000, 13);
    if (csdarange<=0) return false;
    //if (length/csdarange>0.69 && length/csdarange<1.05) return true;
    if (length/csdarange>0.9 && length/csdarange<1.1) return true;
  }

  return false;  
}

//**************************************************
bool MuonPIDACut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  Float_t pida = pdAnaUtils::ComputePIDA(*static_cast<AnaParticlePD*>(box.MainTrack));  
  if (fabs(pida-15.)<4) return true;
  else return false;
}


//**************************************************
void stoppingMuonSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void stoppingMuonSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

