#ifndef LifetimeVariation_h
#define LifetimeVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "CalorimetryAlg.hxx"

class LifetimeVariation: public EventVariationBase, public BinnedParams{
public:
  
  /// Instantiate the PID systematic. nbins is the number of
  /// bins in the PDF
  LifetimeVariation();
  
  virtual ~LifetimeVariation(){} 
  
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
