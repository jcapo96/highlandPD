#ifndef ResidualRangeVariation_h
#define ResidualRangeVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"

class ResidualRangeVariation: public EventVariationBase, public BinnedParams{
public:

  /// This systematic has 0 parameters
  ResidualRangeVariation();
  
  virtual ~ResidualRangeVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:

  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
};

#endif
