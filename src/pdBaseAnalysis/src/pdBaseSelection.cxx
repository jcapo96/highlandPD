#include "pdBaseSelection.hxx"
#include "BasicUtils.hxx"

//********************************************************************
pdBaseSelection::pdBaseSelection(bool forceBreak):SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************  
}

//********************************************************************
void pdBaseSelection::DefineSteps(){
//********************************************************************

  AddStep(StepBase::kCut, "event quality",      new EventQualityCut(),           true);


  SetBranchAlias(0,"dummy");
  
  // No preselection for the moment
  SetPreSelectionAccumLevel(-1);

}

//**************************************************
bool EventQualityCut::Apply(AnaEventC& eventC, ToyBoxB& box) const{
//**************************************************

  (void)box;

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  if(event.GetIsMC())  return true;              // This is MC, ignore DQ. 

  if (enableDQCut) {
    if(!event.DataQuality->GoodDaq ) return false;  // Bad Detector Data quality
  }

  if (enableBeamQualityCut) {
    if(!event.Beam->GoodSpill ) return false;      // Bad Spill
  }

  return true;
}

//**************************************************
void pdBaseSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

}


