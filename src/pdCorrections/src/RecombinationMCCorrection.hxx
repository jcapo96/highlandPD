#ifndef RecombinationMCCorrection_h
#define RecombinationMCCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"
#include "BinnedParams.hxx"

/// This class corrects MC tracks dEdx for different
/// recombination parameters

class RecombinationMCCorrection: public CorrectionBase {

public:
  
  RecombinationMCCorrection();

  virtual ~RecombinationMCCorrection();

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);

protected:

  BinnedParams* _params;
  Calorimetry* _cal;

};

#endif
