#include "pionMomentumSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
pionMomentumSelection::pionMomentumSelection(bool forceBreak)
    : SelectionBase(forceBreak, EventBoxId::kEventBoxPD) {
//********************************************************************
}

//********************************************************************
void pionMomentumSelection::DefineSteps() {
//********************************************************************
  AddStep(StepBase::kCut, "pass all events", new pionMomentumPassAllCut());
  SetBranchAlias(0, "trunk");
}

//********************************************************************
void pionMomentumSelection::InitializeEvent(AnaEventC& eventBB) {
//********************************************************************
  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB);
  if (!event.EventBoxes[EventBoxId::kEventBoxPD]) {
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();
  }
  boxUtils::FillLongTracks(event, static_cast<SubDetId::SubDetEnum>(GetDetectorFV()));
}

//********************************************************************
void pionMomentumSelection::DefineDetectorFV() {
//********************************************************************
  SetDetectorFV(SubDetId::kSubdet1_1);
}

//********************************************************************
bool pionMomentumPassAllCut::Apply(AnaEventC& /*event*/, ToyBoxB& /*box*/) const {
//********************************************************************
  return true;
}
