#include "LifetimeVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>

//********************************************************************
LifetimeVariation::LifetimeVariation(): HitVariationBase(),BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","Lifetime", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void LifetimeVariation::Apply(const ToyExperiment& toy, AnaHitPD& hit){
//********************************************************************

  double fEventT0=500;

  // We need the errors part of the data file but as well the relative uncertainty for sigma
  Float_t mean_corr, mean_var;
  Int_t mean_index;
  
  // Get the systematics parameters
  if (!GetBinValues(0.5, mean_corr,  mean_var,  mean_index))  return;
    
  // Set the nominal electron lifetime and get the lifetime correction
  fCalo.SetElectronLifetime(mean_corr);
  Float_t corr_nominal = fCalo.LifetimeCorrection(hit.PeakTime, fEventT0);
  
  // Set the varied lifetime and get the correction
  fCalo.SetElectronLifetime(mean_corr +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);        
  Float_t corr = fCalo.LifetimeCorrection(hit.PeakTime, fEventT0);
  
  // Uncorrect for the nominal and apply the varied correction only to dQdx
  hit.dQdx      *= corr/corr_nominal;
  hit.dQdx_corr *= corr/corr_nominal;
}

//********************************************************************
bool LifetimeVariation::UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit){
  //********************************************************************
  
  hit.dQdx      = original_hit.dQdx;
  hit.dQdx_corr = original_hit.dQdx_corr;

  // Don't reset the spill to corrected
  return false;
}
