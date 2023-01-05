#ifndef CalorimetryCalibration_h
#define CalorimetryCalibration_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects each hit of a reconstructed track
/// It modifies the pitch by SCE and recomputes dQdx

class CalorimetryCalibration: public CorrectionBase {

public:
  
  CalorimetryCalibration();

  virtual ~CalorimetryCalibration() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  SpaceCharge* _sce;
  Calorimetry* _cal;
};

#endif
