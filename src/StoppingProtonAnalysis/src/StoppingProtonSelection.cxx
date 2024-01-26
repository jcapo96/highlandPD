#include "StoppingProtonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
StoppingProtonSelection::StoppingProtonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void StoppingProtonSelection::DefineSteps(){
//********************************************************************

  //basic cuts → see pdBaseAnalysis/src/pdBaseSelection
  AddStep(StepBase::kAction, "find Pandora track"    , new FindBeamTrackAction(), true);
  AddStep(StepBase::kCut   , "beam track in TPC"     , new BeamTrackExistsCut() , true);
  AddStep(StepBase::kCut   , "beam pdg filter proton", new BeamPDGCut(2212)     , true);
  AddStep(StepBase::kCut   , "BEAM quality cut"      , new BeamQualityCut()     , true);
  
  //proton cuts
  AddStep(StepBase::kCut   , "no daughter cut"       , new NoDaughterCut()      , true);
  AddStep(StepBase::kCut   , "proton CSDA range"     , new CSDARangeCut()       , true);
  
  SetBranchAlias(0,"trunk");
}

//**************************************************
bool NoDaughterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if(!box.MainTrack) return false;
  
  if(box.MainTrack->DaughtersIDs.size() == 0)return true;
  else return false;  
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
void StoppingProtonSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void StoppingProtonSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

