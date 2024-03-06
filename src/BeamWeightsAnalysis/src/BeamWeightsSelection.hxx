#ifndef BeamWeightsSelection_h
#define BeamWeightsSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxdEdx.hxx"
#include "EventBoxdEdx.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class BeamWeightsSelection: public SelectionBase{
 public:
  BeamWeightsSelection(bool forceBreak=true);
  virtual ~BeamWeightsSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxdEdx();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxdEdx::kAllTracksdEdx;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxdEdx::kAllTracksdEdx;return 1;}

  //------------------
  
protected:

};

//---- Steps and actions for the analysis
class BeamPionCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamPionCut();}
};

class NonEmptyEventCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new NonEmptyEventCut();}
};

class EventHasPandoraCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new EventHasPandoraCut();}
};

class DummyCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const {return true;}
  StepBase* MakeClone(){return new DummyCut();}
};

#endif
