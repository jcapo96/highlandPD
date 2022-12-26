#ifndef SCECorrection_h
#define SCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects each hit of a reconstructed track

class SCECorrection: public CorrectionBase {

public:
  
  SCECorrection();

  virtual ~SCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  SpaceCharge* _sce;
  Calorimetry* _cal;
  
};

#endif
