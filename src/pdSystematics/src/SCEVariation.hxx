#ifndef SCEVariation_h
#define SCEVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class SCEVariation: public EventVariationBase, public BinnedParams{
public:

  SCEVariation();
  virtual ~SCEVariation(); 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);

  void VarySCEMapGlobally(const ToyExperiment& toy); //apply a single variation to all the map
  void VarySCEMapLocally(const ToyExperiment& toy);  //apply an independent variation to each voxel

  SpaceCharge* GetToySCE(int itoy) const {return _sce[itoy];}
  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  Calorimetry* _cal;
  SpaceCharge* _sce[100]; //a better way to set the number of toys is clearly needed :)
};

#endif
