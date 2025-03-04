#ifndef secondaryKaonSelectionNoBranches_h
#define secondaryKaonSelectionNoBranches_h

#include "secondaryKaonAnalysis.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxKaon.hxx"
#include "EventBoxKaon.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

//same as the secondary kaon selection but 
//avoiding branches. Useful for dEdx analysis

class secondaryKaonSelectionNoBranches: public SelectionBase{
 public:
  secondaryKaonSelectionNoBranches(bool forceBreak=true);
  virtual ~secondaryKaonSelectionNoBranches(){}

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
  
  static bool IsSuperCandidate(AnaParticlePD* part);
  static bool MuonIsTrack(AnaParticlePD* part);
  static bool MuonChi2(AnaParticlePD* part);
  static bool MuonMom(AnaParticlePD* part);
  static bool KaonMuonAngle(AnaParticlePD* part);
  static bool KaonMuonDistance(AnaParticlePD* part);

protected:

};

class SuperCandidateCut: public StepBase{
 public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new SuperCandidateCut();}
};

#endif
