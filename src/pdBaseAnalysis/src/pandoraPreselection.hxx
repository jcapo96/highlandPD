#ifndef pandoraPreselection_h
#define pandoraPreselection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxPD.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class pandoraPreselection: public SelectionBase{
 public:
  pandoraPreselection(bool forceBreak=true);
  virtual ~pandoraPreselection(){}

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

protected:

  Int_t _FindMainTrackStepIndex;
  Int_t _TotalMultiplicityCutIndex;
};


class FindPandoraTrackAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new FindPandoraTrackAction();}
};

class CandidateExistsCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new CandidateExistsCut();}
};

class BeamQualityCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamQualityCut();}
  bool useIsBeamLike;
  BeamQualityCut(){ useIsBeamLike = (bool)ND::params().GetParameterI("pdBaseAnalysis.UseIsBeamLike");}
};

class CandidateIsBeamCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new CandidateIsBeamCut();}
  bool useIsBeamLike;
  CandidateIsBeamCut(){ useIsBeamLike = (bool)ND::params().GetParameterI("pdBaseAnalysis.UseIsBeamLike");}
};

#endif
