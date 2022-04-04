#ifndef BeamPartIdEffWeight_h
#define BeamPartIdEffWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class BeamPartIdEffWeight: public EventWeightBase, public BinnedParams {
public:
  
  BeamPartIdEffWeight();
  virtual ~BeamPartIdEffWeight(){}

  /// computes the weight for this event
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);

  /// Is this true object relevant for this systematic
  bool IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const;

};

#endif
