#ifndef LifetimeVariation_h
#define LifetimeVariation_h

#include "HitVariationBase.hxx"
#include "BinnedParams.hxx"
#include "CalorimetryAlg.hxx"

class LifetimeVariation: public HitVariationBase, public BinnedParams{
public:
  
  LifetimeVariation();
  
  virtual ~LifetimeVariation(){} 
  
  /// Apply the systematic provided the hit
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit);
  
  /// Undo the variation provided the hit
  bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit);
  
protected:
  calo::CalorimetryAlg fCalo;  
};

#endif
