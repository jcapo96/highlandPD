#ifndef CalorimetryCalibration_h
#define CalorimetryCalibration_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class applies the full calibration chain to
/// MC events

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
