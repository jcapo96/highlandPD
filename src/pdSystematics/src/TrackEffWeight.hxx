#ifndef TrackEffWeight_h
#define TrackEffWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class TrackEffWeight: public EventWeightBase, public BinnedParams {
public:
  
  TrackEffWeight();
  virtual ~TrackEffWeight(){}
  
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);


  bool IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const;

  Int_t GetRelevantTrueObjectGroups(const SelectionBase& sel, Int_t* IDs) const;
  
};

#endif
