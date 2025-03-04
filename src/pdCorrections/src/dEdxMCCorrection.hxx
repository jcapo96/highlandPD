#ifndef dEdxMCCorrection_h
#define dEdxMCCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "BinnedParams.hxx"

#include "TRandom3.h"

/// This class corrects MC dEdx by increasing each MPV and sigma

class dEdxMCCorrection: public CorrectionBase {

public:
  
  dEdxMCCorrection();

  virtual ~dEdxMCCorrection();

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  BinnedParams* _params_1;
  BinnedParams* _params_2;
  BinnedParamsParams* _paramsparams_1;
  BinnedParamsParams* _paramsparams_2;

  TRandom3* _random;

};

#endif
