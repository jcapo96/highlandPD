#include "RecombinationVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>

// MC:   ModBoxA = 0.937±0.020 (ArgoNeuT: 0.93±0.02), ModBoxB = 0.209±0.006 (ArgoNeuT: 0.212±0.002)
// DATA: ModBoxA = 0.854±0.017 (ArgoNeuT: 0.93±0.02) ,ModBoxB = 0.208±0.006 (ArgoNeuT: 0.212±0.002)
                                                         
//********************************************************************
RecombinationVariation::RecombinationVariation(): HitVariationBase(true){
//********************************************************************

  _modBoxA = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","RecombinationA", BinnedParams::k1D_SYMMETRIC);
  _modBoxB = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","RecombinationB", BinnedParams::k1D_SYMMETRIC);
  
  // Read the systematic source parameters from the data files
  SetNParameters(_modBoxA->GetNBins()+_modBoxB->GetNBins());
}

//********************************************************************
void RecombinationVariation::Apply(const ToyExperiment& toy, AnaEventC& eventC){
//********************************************************************

  /* Since the systematic source parameters are hit/particle independent we can set them at the event level */
  
  // We need the errors part of the data file but as well the relative uncertainty for sigma
  Float_t A_mean, A_var, B_mean, B_var;
  Int_t A_index, B_index;
  
  // Get the systematics parameters (This can be outside the loop if it does not depend on Particle properties)
  if (!_modBoxA->GetBinValues(0.5, A_mean,  A_var,  A_index))  return;
  if (!_modBoxB->GetBinValues(0.5, B_mean,  B_var,  B_index))  return;
  
  Float_t ModBoxA = A_mean +  A_var*toy.GetToyVariations(_index)->Variations[A_index];
  Float_t ModBoxB = B_mean +  B_var*toy.GetToyVariations(_index)->Variations[B_index+_modBoxA->GetNBins()]; 
  
  // Set the varied parameters
  fCalo.SetModBoxA(ModBoxA);
  fCalo.SetModBoxB(ModBoxB);        
}

//********************************************************************
void RecombinationVariation::Apply(const ToyExperiment& toy, AnaHitPD& hit){
//********************************************************************
            
  // compute the dEdx from the varied dQdx
  hit.dEdx      = fCalo.dEdx_from_dQdx(hit.dQdx);
  hit.dEdx_corr = fCalo.dEdx_from_dQdx(hit.dQdx_corr);
}

//********************************************************************
bool RecombinationVariation::UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit){
  //********************************************************************
  
  hit.dEdx      = original_hit.dEdx;
  hit.dEdx_corr = original_hit.dEdx_corr;

  // Don't reset the spill to corrected
  return false;
}
