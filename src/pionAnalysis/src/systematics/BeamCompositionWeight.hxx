#ifndef BeamCompositionWeight_h
#define BeamCompositionWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class BeamCompositionWeight: public EventWeightBase, public BinnedParams {
public:
  
  BeamCompositionWeight();
  virtual ~BeamCompositionWeight(){}
  
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);

};

#endif
