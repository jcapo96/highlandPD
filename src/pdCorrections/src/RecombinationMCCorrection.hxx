#ifndef RecombinationMCCorrection_h
#define RecombinationMCCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"
#include "BinnedParams.hxx"

/// This class defines a correction that affects each hit of a reconstructed track
/// It modifies the pitch by SCE and recomputes dQdx

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
