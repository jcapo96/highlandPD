#ifndef ParticlePositionSCECorrection_h
#define ParticlePositionSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

/// This class corrects tracks' information by SCE effect.
/// If tracks have trajectory points, they are corrected
/// individually and then track length and position start/end
/// recomputed. If tracks don't have trajectory points,
/// only position start/end are corrected.

class ParticlePositionSCECorrection: public CorrectionBase {

public:
  
  ParticlePositionSCECorrection();

  virtual ~ParticlePositionSCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  SpaceCharge* _sce;
};

#endif
