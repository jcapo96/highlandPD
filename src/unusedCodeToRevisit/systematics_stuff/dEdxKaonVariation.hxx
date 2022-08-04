#ifndef dEdxKaonVariation_h
#define dEdxKaonVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "pdCalorimetryUtils.hxx"

class dEdxKaonVariation: public EventVariationBase, public BinnedParams{
public:

  /// This systematic has 0 parameters
  dEdxKaonVariation();
  
  virtual ~dEdxKaonVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:

  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
};

#endif
