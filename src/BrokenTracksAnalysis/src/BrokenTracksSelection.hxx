#ifndef BrokenTracksSelection_h
#define BrokenTracksSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxdEdx.hxx"
#include "EventBoxdEdx.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class BrokenTracksSelection: public SelectionBase{
 public:
  BrokenTracksSelection(bool forceBreak=true);
  virtual ~BrokenTracksSelection(){}

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
class BeamMuonCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamMuonCut();}
};

class BeamTrackQualityCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamTrackQualityCut();}
};

class TrackReachesAPABorderCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new TrackReachesAPABorderCut();}
};

class TrackEndsInAPABorderCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new TrackEndsInAPABorderCut();}
};

class TrackIsBrokenCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new TrackIsBrokenCut();}
};

#endif
