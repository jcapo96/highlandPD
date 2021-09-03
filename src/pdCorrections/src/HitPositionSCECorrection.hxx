#ifndef HitPositionSCECorrection_h
#define HitPositionSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "pdSpaceCharge.hxx"

/// This class defines a correction that affects each hit of a reconstructed track

class HitPositionSCECorrection: public CorrectionBase {

public:
  
  HitPositionSCECorrection();

  virtual ~HitPositionSCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  pdspacecharge::pdSpaceCharge sce;
  
};

#endif
