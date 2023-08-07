#ifndef KaonEnergyCorrection_h
#define KaonEnergyCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "BinnedParams.hxx"

/// This class defines a correction that affects end/start position of each track

class KaonEnergyCorrection: public CorrectionBase {

public:
  
  KaonEnergyCorrection();

  virtual ~KaonEnergyCorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  BinnedParams* _params;

};

#endif
