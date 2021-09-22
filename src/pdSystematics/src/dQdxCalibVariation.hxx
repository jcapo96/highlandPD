#ifndef dQdxCalibVariation_h
#define dQdxCalibVariation_h

#include "HitVariationBase.hxx"
#include "BinnedParams.hxx"

class dQdxCalibVariation: public HitVariationBase, public BinnedParams{
public:
  
  dQdxCalibVariation();
  
  virtual ~dQdxCalibVariation(){} 

  /// Apply the systematic provided the hit
  using HitVariationBase::Apply;
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit);
  
  /// Undo the variation provided the hit
  using HitVariationBase::UndoSystematic;
  bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit);

protected:
  BinnedParams* _paramX;
  BinnedParams* _paramYZ;
  BinnedParams* _paramNQ;
  BinnedParams* _paramCcal;
};

#endif
