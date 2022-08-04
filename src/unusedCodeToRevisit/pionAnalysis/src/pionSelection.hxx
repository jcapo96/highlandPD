#ifndef pionSelection_h
#define pionSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxPD.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"
#include "CNNUtils.hxx"
#include "Parameters.hxx"

/* 
   This file contains the Pion Selection developed b J. Calcutt and F. Stocker, relized in the class pionSelection, which is a collection of steps (StepBase class). 
   The steps (cuts or actions) are added to the pionSelection class in the source file. 
*/


class pionSelection: public SelectionBase{
 public:
  pionSelection(bool forceBreak=true);
  virtual ~pionSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxPD();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxPD::kCandidateAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxPD::kCandidateAndDaughters;return 1;}

  //------------------

protected:

  // TODO. Not used
  Int_t _PionRangeCutIndex;
  Int_t _PionRangeStepIndex;
};


class BeamPionCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamPionCut();}
};

class CandidateIsTrackCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new CandidateIsTrackCut();}
};

class PionEndsAPA3Cut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new PionEndsAPA3Cut();}
};

class PionPassChi2Cut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new PionPassChi2Cut();}
};

class NoPionDaughterCut: public StepBase{
 public:
  NoPionDaughterCut(){
    _cnnRecomputeCut = ND::params().GetParameterD("pionAnalysis.Systematics.CNNRecomputeCut");
    _cnnSystematicEnabled = ND::params().GetParameterD("pionAnalysis.Systematics.EnableCNN");
  }
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new NoPionDaughterCut();}
  Float_t _cnnRecomputeCut;
  bool _cnnSystematicEnabled;
};

class PionCexCut: public StepBase{
 public:
  PionCexCut(){
    _cnnRecomputeCut = ND::params().GetParameterD("pionAnalysis.Systematics.CNNRecomputeCut");
    _cnnSystematicEnabled = ND::params().GetParameterD("pionAnalysis.Systematics.EnableCNN");
  }
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new PionCexCut();}
  Float_t _cnnRecomputeCut;
  bool _cnnSystematicEnabled;
};

class PionAbsorptionCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new PionAbsorptionCut();}
  PionCexCut _cut;
};




#endif
