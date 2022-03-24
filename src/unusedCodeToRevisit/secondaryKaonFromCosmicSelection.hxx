#ifndef secondaryKaonFromCosmicSelection_h
#define secondaryKaonFromCosmicSelection_h

#include "secondaryKaonAnalysis.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxKaon.hxx"
#include "EventBoxPD.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"


class secondaryKaonFromCosmicSelection: public SelectionBase{
 public:
  secondaryKaonFromCosmicSelection(bool forceBreak=true);
  virtual ~secondaryKaonFromCosmicSelection(){}

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
class GetKaonsFromCosmicsAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetKaonsFromCosmicsAction();}
};

#endif
