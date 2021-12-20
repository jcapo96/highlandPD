#ifndef secondaryKaonSelection_h
#define secondaryKaonSelection_h

#include "secondaryKaonAnalysis.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxKaon.hxx"
#include "EventBoxPD.hxx"
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
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxPD::kCandidateAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxPD::kCandidateAndDaughters;return 1;}

  //------------------

protected:

};

//---- Steps and actions for the analysis
class BeamFilterCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamFilterCut();}
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
  double _cut_value;
  MuonChi2Cut(double cut_value){_cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonChi2Cut(_cut_value);}
};

class MuonCNNCut: public StepBase{
 public:
  double _cut_value;
  MuonCNNCut(double cut_value){_cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonCNNCut(_cut_value);}
};

class MuonRangeMomCut: public StepBase{
 public:
  double _cut_value;
  double _mean_value;
  MuonRangeMomCut(double mean_value,double cut_value){_mean_value = mean_value; _cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuonRangeMomCut(_mean_value,_cut_value);}
};

class KaonCNNCut: public StepBase{
 public:
  double _cut_value;
  KaonCNNCut(double cut_value){_cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonCNNCut(_cut_value);}
};

class KaonMuonAngleCut: public StepBase{
 public:
  double _cut_value;
  KaonMuonAngleCut(double cut_value){_cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonMuonAngleCut(_cut_value);}
};

class KaonMuonDistanceCut: public StepBase{
 public:
  double _cut_value;
  KaonMuonDistanceCut(double cut_value){_cut_value = cut_value;};
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new KaonMuonDistanceCut(_cut_value);}
};

#endif
