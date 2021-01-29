#ifndef dQdxCalibVariation_h
#define dQdxCalibVariation_h

#include "HitVariationBase.hxx"
#include "BinnedParams.hxx"

class dQdxCalibVariation: public HitVariationBase, public BinnedParams{
public:
  
  dQdxCalibVariation();
  
  virtual ~dQdxCalibVariation(){} 

  /// Apply the systematic provided the hit
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit);
  
  /// Undo the variation provided the hit
  bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit);

protected:
  BinnedParams* _calibX;
  BinnedParams* _calibYZ;
};

#endif
