#ifndef RecombinationVariation_h
#define RecombinationVariation_h

#include "HitVariationBase.hxx"
#include "BinnedParams.hxx"
#include "pdCalorimetryUtils.hxx"

class RecombinationVariation: public HitVariationBase, public BinnedParams{
public:
  
  RecombinationVariation();
  
  virtual ~RecombinationVariation(){} 

  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Apply the systematic provided the hit
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit);
  
  /// Undo the variation provided the hit
  using HitVariationBase::UndoSystematic;
  bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit);

protected:
  pdCalorimetryUtils fCalo;

  BinnedParams* _modBoxA;
  BinnedParams* _modBoxB;  
};

#endif
