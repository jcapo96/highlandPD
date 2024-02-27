#include "BeamWeightsSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
BeamWeightsSelection::BeamWeightsSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxdEdx) {
//********************************************************************

}

//********************************************************************
void BeamWeightsSelection::DefineSteps(){
//********************************************************************

  AddStep(StepBase::kCut,    "dummy cut",           new DummyCut(), true);

  //split the selection in branches, one for kaons, one for protons and another for pions/muons/electrons
  AddSplit(3);

  //pions/muons/electrons
  AddStep(0, StepBase::kCut, "beam pions",        new BeamPionCut(),      true);
  AddStep(0, StepBase::kCut, "event not empty",   new NonEmptyEventCut(), true);
  AddStep(0, StepBase::kCut, "event has pandora", new EventHasPandoraCut(), true);

  //protons
  AddStep(1, StepBase::kCut, "beam protons",      new BeamPDGCut(2212),   true);
  AddStep(1, StepBase::kCut, "event not empty",   new NonEmptyEventCut(), true);
  AddStep(1, StepBase::kCut, "event has pandora", new EventHasPandoraCut(), true);
  
  //kaons
  AddStep(2, StepBase::kCut, "beam kaons",        new BeamPDGCut(321),    true);
  AddStep(2, StepBase::kCut, "event not empty",   new NonEmptyEventCut(), true);
  AddStep(2, StepBase::kCut, "event has pandora", new EventHasPandoraCut(), true);

  SetBranchAlias(0,"beam pions",0);
  SetBranchAlias(1,"beam protons",1);
  SetBranchAlias(2,"beam kaons",2);
}

//**************************************************
bool BeamPionCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut beamCut0(11);
  BeamPDGCut beamCut1(13);
  BeamPDGCut beamCut2(211);

  if(beamCut0.Apply(event,boxB) ||
     beamCut1.Apply(event,boxB) || beamCut2.Apply(event,boxB))
    return true;
  else return false;
}

//**************************************************
bool NonEmptyEventCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(static_cast<AnaEvent*>(&event)->EventInfo);

  return !info->EmptyEvent;
}

//**************************************************
bool EventHasPandoraCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(static_cast<AnaEvent*>(&event)->EventInfo);

  return info->HasPandora;
}

//**************************************************
void BeamWeightsSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();
  
  boxUtils::FillCandidateAndDaughters(event);
}
