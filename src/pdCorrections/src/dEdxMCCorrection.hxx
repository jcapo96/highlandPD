#ifndef dEdxMCCorrection_h
#define dEdxMCCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "BinnedParams.hxx"

/// This class defines a correction that affects end/start position of each track

class dEdxMCCorrection: public CorrectionBase {

public:
  
  dEdxMCCorrection();

  virtual ~dEdxMCCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  BinnedParams* _params;

};

#endif
