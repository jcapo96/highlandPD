#include "stoppingProtonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"
#include "pdSystId.hxx"

//********************************************************************
stoppingProtonSelection::stoppingProtonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void stoppingProtonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")
  //AddStep(StepBase::kAction, "find main track",     new FindBeamTrackAction());  
  AddStep(StepBase::kAction, "find Pandora track",  new FindBeamTrackAction());
  AddStep(StepBase::kCut,    "beam track in TPC",   new BeamPDGCut(2212));
  //AddStep(StepBase::kCut,    "pandora reco worked", new CandidateIsBeamCut());
  AddStep(StepBase::kCut,    "proton CSDA range",   new ProtonCSDARangeCut());
  
  SetBranchAlias(0,"trunk");
}

//**************************************************
bool BeamProtonAngleCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  //Cast the particle and define vector
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.MainTrack);
  double beamdir[3] = {-0.18,-0.20,0.96};//beam direction
  double partdir[3] = {part->DirectionStart[0],part->DirectionStart[1],part->DirectionStart[2]};

  //Cut
  if((beamdir[0]*partdir[0]+beamdir[1]*partdir[1]+beamdir[2]*partdir[2])>0.93)return true;
  else return false;  
}

//**************************************************
bool ProtonCSDARangeCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  //check if it exists a beam particle
  AnaBeamPD* beam         = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;
  if(beamPart){
    Float_t mom = beamPart->Momentum;
    //if (beamPart->TrueObject) mom =  static_cast<AnaTrueParticle*>(beamPart->TrueObject)->Momentum;
    Float_t length = static_cast<AnaParticlePD*>(box.MainTrack)->Length;
    Float_t csdarange = pdAnaUtils::ComputeCSDARange(mom*1000, 2212);
    if (csdarange<=0) return false;
    //if (length/csdarange>0.69 && length/csdarange<1.05) return true;
    if (length/csdarange>0.74 && length/csdarange<1.09) return true;
  }

  return false;  
}

//**************************************************
bool ProtonPIDACut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
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
void stoppingProtonSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void stoppingProtonSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

//**************************************************
bool stoppingProtonSelection::IsRelevantTrueObjectForSystematicInToy(const AnaEventC& event, const ToyBoxB& boxB, AnaTrueObjectC* trueObj, SystId_h systId, Int_t branch) const{
//**************************************************

  (void)event;
  (void)branch;

  // Cast the ToyBox to the appropriate type
  const ToyBoxPD& box = *static_cast<const ToyBoxPD*>(&boxB); 


  AnaTrueParticleB* trueTrack = static_cast<AnaTrueParticleB*>(trueObj);
  
  
  // Apply to all objects,  not fine-tuning
  //  if (!systTuning::APPLY_SYST_FINE_TUNING) return true;
 
  // Main track "mode",  will only consider certain true tracks of interest

  if(systId == pdSystId::kTrackEff){
    if (box.MainTrack->GetTrueParticle()){
      // At first order the inclusive selection only depends on the tracking efficiency of the proton candidate.
      if (trueTrack->ID  == box.MainTrack->GetTrueParticle()->ID) return true; 
    
      // Consider also the case in which the muon candidate is not a true muon but this track it is
      if (trueTrack->PDG == 2212 && box.MainTrack->GetTrueParticle()->PDG!=13) return true;     
      
      return false;

    } 
    // Apply for muon 
    else if (trueTrack->PDG == 2212)  
      return true;
    
    return false;
  } 

  return true;
}
