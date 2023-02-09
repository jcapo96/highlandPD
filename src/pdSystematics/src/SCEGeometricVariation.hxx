#ifndef SCEGeometricVariation_h
#define SCEGeometricVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class SCEGeometricVariation: public EventVariationBase, public BinnedParams{
public:

  SCEGeometricVariation();
  virtual ~SCEGeometricVariation(); 

  void Initialize();
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);

  void VarySCEMap(const ToyExperiment& toy);
  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  SpaceCharge* _sce[100];
};

#endif
