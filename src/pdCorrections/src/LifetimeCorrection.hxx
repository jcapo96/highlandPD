#ifndef LifetimeCorrection_h
#define LifetimeCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects each hit of a reconstructed track

class LifetimeCorrection: public CorrectionBase {

public:
  
  LifetimeCorrection();

  virtual ~LifetimeCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  //this probably needs a reimplementation, it cannot be initialized here. Possibly this should be a singleton.
  Calorimetry* _cal;
  
};

#endif
