#ifndef HitPitchSCECorrection_h
#define HitPitchSCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"

/// This class defines a correction that affects each hit of a reconstructed track

class HitPitchSCECorrection: public CorrectionBase {

public:
  
  HitPitchSCECorrection();

  virtual ~HitPitchSCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  //this probably needs a reimplementation, it cannot be initialized here. Possibly this should be a singleton.
  //edit: it isn't really a problem to define here a spacecharge instance that is initialized afterwards in the
  //'initialize' method. The problem will come when it needs to be used in more places than here. Do we want
  //to have an instance for each place in which it has to be used, or a single instance for all of them? Second
  //option would be better if there are effects across the usages. To be studied.
  SpaceCharge* _sce;
  
};

#endif
