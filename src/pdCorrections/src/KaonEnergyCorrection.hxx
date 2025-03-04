#ifndef KaonEnergyCorrection_h
#define KaonEnergyCorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "BinnedParams.hxx"

/// This class modifies kaons dEdx
/// deprecated, to be removed

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
