#ifndef dEdxUtils_h
#define dEdxUtils_h

#include "TF1.h"

#include "dQdxYZCalibration.hxx"
#include "SpaceCharge.hxx"

namespace dEdxUtils{

  Double_t Langaus(Double_t *x, Double_t *par);

  void SetFitParameters(TF1* f, Double_t max);

  Double_t Recombination(Double_t E);

  Double_t ElectricField(const SpaceCharge* sce, Double_t x, Double_t y, Double_t z);

  bool IsInterestingHit(AnaHitPD& hit);
}

#endif
