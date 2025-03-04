#ifndef DaughterFinderCorrection_h
#define DaughterFinderCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// Not sure about this, seems secondary kaon XS stuff

class DaughterFinderCorrection: public CorrectionBase {

public:
  
  DaughterFinderCorrection();

  virtual ~DaughterFinderCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

};

#endif
