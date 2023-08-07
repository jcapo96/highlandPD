#ifndef ProtonBackgroundWeight_h
#define ProtonBackgroundWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class ProtonBackgroundWeight: public EventWeightBase, public BinnedParams {
public:
  
  ProtonBackgroundWeight();
  virtual ~ProtonBackgroundWeight(){}
  
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);


  /// Check whether a AnaRecObject is relevant for this systematic or not
  virtual bool IsRelevantRecObject(const AnaEventC&, const AnaRecObjectC&) const;
  
};

#endif
