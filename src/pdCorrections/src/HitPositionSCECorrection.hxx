#ifndef HitPositionSCECorrection_h
#define HitPositionSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"

/// This class defines a correction that affects each hit of a reconstructed track

class HitPositionSCECorrection: public CorrectionBase {

public:
  
  HitPositionSCECorrection();

  virtual ~HitPositionSCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  //this probably needs a reimplementation, it cannot be initialized here. Possibly this should be a singleton.
  SpaceCharge* sce;
  
};

#endif
