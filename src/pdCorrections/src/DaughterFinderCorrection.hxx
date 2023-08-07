#ifndef DaughterFinderCorrection_h
#define DaughterFinderCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "Calorimetry.hxx"

/// This class defines a correction that affects each hit of a reconstructed track
/// It modifies the pitch by SCE and recomputes dQdx

class DaughterFinderCorrection: public CorrectionBase {

public:
  
  DaughterFinderCorrection();

  virtual ~DaughterFinderCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

};

#endif
