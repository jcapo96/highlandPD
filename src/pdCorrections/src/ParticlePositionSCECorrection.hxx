#ifndef ParticlePositionSCECorrection_h
#define ParticlePositionSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects end/start position of each track

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
