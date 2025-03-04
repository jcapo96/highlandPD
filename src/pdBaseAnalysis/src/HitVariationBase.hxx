#ifndef HitVariationBase_h
#define HitVariationBase_h

#include "EventVariationBase.hxx"
#include "pdDataClasses.hxx"

///Migue: Anselmo's stuff, probably to be removed

/// This is the base class for systematics applied at the hit level
/// - Variation: modify some aspect of the input data. The selection needs to be redone 
///   in general. This is what we call Toy Experiments

class HitVariationBase: public EventVariationBase{
  
public:
  
  /// Create the systematic, specifying the number of systematic parameters and whether dEdx is recomputed inside this variation
  HitVariationBase(bool dEdx=false, UInt_t npar=0):EventVariationBase(npar){dEdx_recomputed=dEdx;}
  
  /// Everyone should have a destructor.
  virtual ~HitVariationBase() {}

  /// By default no neew to do anyting at the event level
  /// This CAN be overridden in the derived class.
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event){}
  
  //----------------------- MANDATORY METHODS -----------------------------------------

  /// Apply the systematic
  /// This MUST be overridden in the derived class.
  virtual void Apply(const ToyExperiment& toy, AnaHitPD& hit) = 0;

  //----------------------- OPTIONAL METHODS -----------------------------------------

  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  /// This undos the variations. If it return true the Spill will be reset  
  /// By default the Spill is reset
  virtual bool UndoSystematic(AnaEventC&){return false;}
  
  /// Undo the variation provided the hit
  virtual bool UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit){return false;}

public:

  /// Indicates whether dEdx is recomputed inside this variation
  bool dEdx_recomputed;

};

#endif
