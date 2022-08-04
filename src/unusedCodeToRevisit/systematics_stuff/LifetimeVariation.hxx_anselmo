#ifndef LifetimeVariation_h
#define LifetimeVariation_h

#include "HitVariationBase.hxx"
#include "BinnedParams.hxx"
#include "pdCalorimetryUtils.hxx"

class LifetimeVariation: public HitVariationBase, public BinnedParams{
public:
  
  LifetimeVariation();
  
  virtual ~LifetimeVariation(){} 
  
  /// Apply the systematic provided the hit
  using HitVariationBase::Apply;
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit);
  
  /// Undo the variation provided the hit
  using HitVariationBase::UndoSystematic;
  bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit);
  
protected:
  pdCalorimetryUtils fCalo;  
};

#endif
