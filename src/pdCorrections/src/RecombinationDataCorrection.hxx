#ifndef RecombinationDataCorrection_h
#define RecombinationDataCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class corrects data tracks dEdx for different
/// recombination parameters

class RecombinationDataCorrection: public CorrectionBase {

public:
  
  RecombinationDataCorrection();

  virtual ~RecombinationDataCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

};

#endif
