#ifndef neutralKaonSelection_h
#define neutralKaonSelection_h

#include "neutralKaonAnalysis.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxNeutralKaon.hxx"
#include "EventBoxKaon.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class neutralKaonSelection: public SelectionBase{
 public:
 neutralKaonSelection(bool forceBreak=true);
  virtual ~neutralKaonSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(); //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxNeutralKaon();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough in many cases
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  //Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}
  //Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}

  //------------------

protected:

};

class GetNeutralKaonsAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new GetNeutralKaonsAction();}
  };

class EventHasVertex: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new EventHasVertex();}
  };

// New Action classes for Preliminary K0 Selection
class FindBeamDaughtersAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new FindBeamDaughtersAction();}
  };

class CheckK0InTruthAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new CheckK0InTruthAction();}
  };

// New Cut classes for Preliminary K0 Selection
class HasEnoughBeamDaughtersCut: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new HasEnoughBeamDaughtersCut();}
  };

class BeamDaughtersDistanceCut: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new BeamDaughtersDistanceCut(_lower_cut, _upper_cut);}
  };

class BeamDaughtersPairDistanceCut: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new BeamDaughtersPairDistanceCut(_lower_cut, _upper_cut);}
  };

#endif
