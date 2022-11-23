#ifndef secondaryKaonSelection_h
#define secondaryKaonSelection_h

#include "secondaryKaonAnalysis.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxKaon.hxx"
#include "EventBoxKaon.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class secondaryKaonSelection: public SelectionBase{
 public:
  secondaryKaonSelection(bool forceBreak=true);
  virtual ~secondaryKaonSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxKaon();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}

  //------------------
  
protected:

};

//---- Steps and actions for the analysis
class BeamHadronCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamHadronCut();}
};

class GetKaonsAction: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetKaonsAction();}
};

class EventHasKaonCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new EventHasKaonCut();}
};

class MuonIsTrackCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonIsTrackCut();}
};

class MuonChi2Cut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonChi2Cut();}
};

class MuonCNNCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonCNNCut();}
};

class MuonRangeMomCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonRangeMomCut();}
};

class KaonCNNCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonCNNCut();}
};

class KaonMuonAngleCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonMuonAngleCut();}
};

class KaonMuonDistanceCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonMuonDistanceCut();}
};

class ProtonChi2Cut: public StepBase{
 public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new ProtonChi2Cut();}
};

#endif
