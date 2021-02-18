#ifndef CNNVariation_h
#define CNNVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "CalorimetryAlg.hxx"

class CNNVariation: public EventVariationBase, BinnedParams{
public:

  /// This systematic has 0 parameters
  CNNVariation();
  
  virtual ~CNNVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:

  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;

  // needed to recompute dEdx
  calo::CalorimetryAlg fCalo;

  Float_t _cnnRecomputeCut;
  
};

#endif
