#ifndef TrackEffWeight_h
#define TrackEffWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class TrackEffWeight: public EventWeightBase, public BinnedParams {
public:
  
  TrackEffWeight();
  virtual ~TrackEffWeight(){}

  /// computes the weight for this event
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);

  /// Is this true object relevant for this systematic
  bool IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const;

};

#endif
