#ifndef HitPitchSCECorrection_h
#define HitPitchSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects each hit of a reconstructed track
/// It modifies the pitch by SCE and recomputes dQdx

class HitPitchSCECorrection: public CorrectionBase {

public:
  
  HitPitchSCECorrection();

  virtual ~HitPitchSCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  Calorimetry* _cal;
};

#endif
