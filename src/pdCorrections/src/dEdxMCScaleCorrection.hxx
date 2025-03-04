#ifndef dEdxMCScaleCorrection_h
#define dEdxMCScaleCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "BinnedParams.hxx"

/// This class defines a correction that scales
/// all MC dEdx values

class dEdxMCScaleCorrection: public CorrectionBase {

public:
  
  dEdxMCScaleCorrection();

  virtual ~dEdxMCScaleCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  BinnedParams* _params;

};

#endif
