#ifndef LifetimeVariation_h
#define LifetimeVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class LifetimeVariation: public EventVariationBase, public BinnedParams{
public:

  LifetimeVariation();
  virtual ~LifetimeVariation();
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  Calorimetry* _cal;
};

#endif
