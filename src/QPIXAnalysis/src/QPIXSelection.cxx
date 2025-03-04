#include "QPIXSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
QPIXSelection::QPIXSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void QPIXSelection::DefineSteps(){
//********************************************************************

  //basic cuts → see pdBaseAnalysis/src/pdBaseSelection
  AddStep(StepBase::kCut, "Dummy", new DummyCut(), true);
    
  SetBranchAlias(0,"trunk");
}

//**************************************************
bool DummyCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  return true;
}

//**************************************************
void QPIXSelection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillLongTracks(event,static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void QPIXSelection::DefineDetectorFV(){
//********************************************************************

    // The detector in which the selection is applied
    SetDetectorFV(SubDetId::kSubdet1_1);
}

