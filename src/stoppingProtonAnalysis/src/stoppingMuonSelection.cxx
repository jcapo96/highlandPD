#include "stoppingMuonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

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
  AddStep(StepBase::kAction, "find Pandora track",  new FindPandoraTrackAction());
  AddStep(StepBase::kCut,    "beam track in TPC",   new CandidateExistsCut());
  AddStep(StepBase::kCut,    "beam muon",           new BeamMuonCut());
  AddStep(StepBase::kCut,    "pandora reco worked", new CandidateIsBeamCut());
  AddStep(StepBase::kCut,    "muon CSDA range",     new MuonCSDARangeCut());
  
  SetBranchAlias(0,"trunk");
  //AddStep(StepBase::kAction, "find vertex",          new FindVertexAction());  // action from duneExampleAnalysis package
  //dunno why this action was being used. Commented for the moment
}

//**************************************************
bool BeamMuonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)boxB;
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  if (trueBeamPart){
    // for MC
    if (trueBeamPart->PDG==-13 || trueBeamPart->PDG==211) return true;
  }
  else{
    // for real DATA
    for(int i = 0; i < (int)beam->PDGs.size(); i++){
      if (beam->PDGs[i] == 13 || beam->PDGs[i] == 211) return true;
    }
  }
  return false;
  /*(void)boxB;
  
  AnaBeam* beam = static_cast<AnaBeam*>(static_cast<AnaEventB*>(&event)->Beam);

  if (event.GetIsMC()){
    if (beam->BeamParticle){
      if (beam->BeamParticle->TrueObject){       
        if (abs(static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject)->PDG)==13 || abs(static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject)->PDG)==211) return true;
        // If electron veto is off consier also electrons
        if (!_electronVeto && abs(static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject)->PDG)==11) return true;
      }
    }
    return false;
  }
  else{
    if (beam->TOF>145 && beam->TOF<170) return true;
    else return false;    
    }*/
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
    Float_t length = static_cast<AnaParticlePD*>(box.MainTrack)->corrected_Length;
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

