#ifndef NominalBeamMomWeight_h
#define NominalBeamMomWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class NominalBeamMomWeight: public EventWeightBase, public BinnedParams {
public:
  
  NominalBeamMomWeight();
  virtual ~NominalBeamMomWeight(){}
  
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);

};

#endif
