#ifndef pdBaseSelection_h
#define pdBaseSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxPD.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

/// This file defines steps that are common to all analyses.
/// Currently the only cut is one on event quality.


class pdBaseSelection: public SelectionBase{
public:
  pdBaseSelection(bool forceBreak=true);
  virtual ~pdBaseSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxPD();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxPD::kCandidateAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t*, Int_t) const{return 0;}

  //------------------
  
};

/// A cut on event quality. Requires good beam and detector data quality flags.
class EventQualityCut: public StepBase {
public:
  using StepBase::Apply;
  
  EventQualityCut(){
    enableDQCut = (bool) ND::params().GetParameterI("pdBaseAnalysis.EnableDataQualityCut");
    enableBeamQualityCut = (bool) ND::params().GetParameterI("pdBaseAnalysis.EnableBeamQualityCut");
  }
  
  bool enableDQCut;
  bool enableBeamQualityCut;
  
  /// Apply the event quality cut. See EventQualityCut class documentation
  /// for details.
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new EventQualityCut();}
};

class FindBeamTrackAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new FindBeamTrackAction();}
};

class BeamTrackExistsCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamTrackExistsCut();}
};

class BeamPDGCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamPDGCut();}
};

class BeamQualityCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamQualityCut();}
};

#endif
