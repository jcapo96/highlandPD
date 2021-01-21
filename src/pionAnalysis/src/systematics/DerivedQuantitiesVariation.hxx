#ifndef DerivedQuantitiesVariation_h
#define DerivedQuantitiesVariation_h

#include "EventVariationBase.hxx"
#include "CalorimetryAlg.hxx"

class DerivedQuantitiesVariation: public EventVariationBase{
public:
  
  /// Constructor
  DerivedQuantitiesVariation();
  
  virtual ~DerivedQuantitiesVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  

protected:

  /// Is this track relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const;

  /// Get the TrackGroup IDs array for this systematic
  Int_t GetRelevantRecObjectGroups(const SelectionBase& sel,     Int_t* IDs) const; 
  
  calo::CalorimetryAlg fCalo;
};

#endif
