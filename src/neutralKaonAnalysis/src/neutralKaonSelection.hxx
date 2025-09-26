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
class FindAllParticlesAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new FindAllParticlesAction();}
  };

class CheckK0InTruthAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new CheckK0InTruthAction();}
  };


class FindNeutralCandidatesAction: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new FindNeutralCandidatesAction();}
  };

// New Cut classes for Preliminary K0 Selection
class HasEnoughParticlesCut: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new HasEnoughParticlesCut();}
  };


class HasNeutralCandidatesCut: public StepBase{
  public:
    using StepBase::Apply;
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new HasNeutralCandidatesCut();}
  };

class VertexParentCountCut: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    VertexParentCountCut(int min_parents = 1);
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new VertexParentCountCut(_min_parents);}
  private:
    int _min_parents;
  };

class VertexDaughterCountCut: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    VertexDaughterCountCut(int min_daughters = 2);
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new VertexDaughterCountCut(_min_daughters);}
  private:
    int _min_daughters;
  };

class VtxDaughtersAreDaughtersOfVtxParentCut: public StepBase{
  public:
    using StepBase::StepBase;
    using StepBase::Apply;
    VtxDaughtersAreDaughtersOfVtxParentCut();
    bool Apply(AnaEventC& event, ToyBoxB& box) const;
    StepBase* MakeClone(){return new VtxDaughtersAreDaughtersOfVtxParentCut();}
  };
#endif
