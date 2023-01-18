#ifndef ParticlePositionSCE_h
#define ParticlePositionSCE_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects end/start position of each track

class ParticlePositionSCE: public CorrectionBase {

public:
  
  ParticlePositionSCE();

  virtual ~ParticlePositionSCE() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  SpaceCharge* _sce;
};

#endif
