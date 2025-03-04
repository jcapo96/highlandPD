#ifndef dEdxVariation_h
#define dEdxVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class dEdxVariation: public EventVariationBase, public BinnedParams{
public:

  dEdxVariation();
  virtual ~dEdxVariation(){} 
  
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
