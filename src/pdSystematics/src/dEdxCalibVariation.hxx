#ifndef dEdxCalibVariation_h
#define dEdxCalibVariation_h

#include "EventVariationBase.hxx"
#include "pdCalorimetryUtils.hxx"

class dEdxCalibVariation: public EventVariationBase{
public:

  /// This systematic has 0 parameters
  dEdxCalibVariation():EventVariationBase(0){}
  
  virtual ~dEdxCalibVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:

  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;

  // needed to recompute dEdx
  pdCalorimetryUtils fCalo;
};

#endif
