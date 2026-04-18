#ifndef pionMomentumSelection_h
#define pionMomentumSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "EventBoxPD.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class pionMomentumSelection : public SelectionBase {
 public:
  explicit pionMomentumSelection(bool forceBreak = true);
  virtual ~pionMomentumSelection() {}

  void DefineSteps();
  void DefineDetectorFV();
  ToyBoxB* MakeToyBox() { return new ToyBoxPD(); }
  void InitializeEvent(AnaEventC&);

  bool FillEventSummary(AnaEventC&, Int_t*) { return false; }
  SampleId_h GetSampleId() { return UNASSIGNEDID; }
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const {
    IDs[0] = EventBoxPD::kLongTracks;
    return 1;
  }
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t*, Int_t) const { return 0; }

 protected:
};

class pionMomentumPassAllCut : public StepBase {
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone() { return new pionMomentumPassAllCut(); }
};

#endif
