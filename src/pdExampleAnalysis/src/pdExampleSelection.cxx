#include "pdExampleSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
pdExampleSelection::pdExampleSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void pdExampleSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")
  AddStep(StepBase::kAction, "find Pandora track",           new FindBeamTrackAction());
  AddStep(StepBase::kCut,    "beam track in TPC",            new BeamTrackExistsCut());
  //AddStep(StepBase::kCut,    "beam particle is proton-like", new BeamPDGCut(2212);
  //AddStep(StepBase::kCut,    "pandora reco worked",          new BeamQualityCut());
  AddStep(StepBase::kCut,    "proton CSDA range",            new CSDARangeCut());
  
  SetBranchAlias(0,"trunk");
}

//**************************************************
bool CSDARangeCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  
  //check if it exists a beam particle
  AnaBeamPD* beam           = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;
  if(beamPart){
    Float_t mom       = beamPart->Momentum;
    Float_t length    = static_cast<AnaParticlePD*>(box.MainTrack)->Length;
    Float_t csdarange = pdAnaUtils::ComputeCSDARange(mom*1000, 2212);
    if (csdarange<=0) return false;
    if (length/csdarange>0.74 && length/csdarange<1.09) return true;
  }

  return false;  
}


//**************************************************
void pdExampleSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void pdExampleSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

