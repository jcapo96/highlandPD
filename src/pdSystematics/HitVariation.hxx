#ifndef HitVariation_h
#define HitVariation_h

#include "pdSpaceCharge.hxx"
#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"

/// This systematic smears the CT of each TPC track segment
class HitVariation: public EventVariationBase, public BinnedParams{
public:
  
  /// Instantiate the PID systematic. nbins is the number of
  /// bins in the PDF
  HitVariation();
  
  virtual ~HitVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster than reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  

protected:

  /// Is this track relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const;

  /// Get the TrackGroup IDs array for this systematic
  Int_t GetRelevantRecObjectGroups(const SelectionBase& sel,     Int_t* IDs) const; 
  
  //for sce corrections
  pdspacecharge::pdSpaceCharge* sce = new pdspacecharge::pdSpaceCharge();
};

#endif
