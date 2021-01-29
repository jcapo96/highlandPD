#include "dQdxCalibVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
dQdxCalibVariation::dQdxCalibVariation(): HitVariationBase(){
//********************************************************************

  _calibX  = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCalibX",  BinnedParams::k1D_SYMMETRIC);
  _calibYZ = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCalibYZ", BinnedParams::k1D_SYMMETRIC);
  
  // Read the systematic source parameters from the data files
  SetNParameters(_calibX->GetNBins()+_calibYZ->GetNBins());
}

//********************************************************************
void dQdxCalibVariation::Apply(const ToyExperiment& toy, AnaHitPD& hit){
//********************************************************************
    
  // We need the errors part of the data file but as well the relative uncertainty for sigma
  Float_t X_mean, X_var, YZ_mean, YZ_var;
  Int_t X_index, YZ_index;
  
  // Get the systematics parameters
  if ( !_calibX->GetBinValues(0.5, X_mean,  X_var,  X_index))  return;
  if (!_calibYZ->GetBinValues(0.5, YZ_mean, YZ_var, YZ_index))  return;
  
  Float_t C_X  =  X_mean +   X_var*toy.GetToyVariations(_index)->Variations[X_index];
  Float_t C_YZ = YZ_mean +  YZ_var*toy.GetToyVariations(_index)->Variations[X_index+_calibX->GetNBins()]; 
    
  hit.dQdx       *= C_X*C_YZ;
  hit.dQdx_NoSCE *= C_X*C_YZ;
}

//********************************************************************
bool dQdxCalibVariation::UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit){
  //********************************************************************
  
  hit.dQdx       = original_hit.dQdx;
  hit.dQdx_NoSCE = original_hit.dQdx_NoSCE;

  // Don't reset the spill to corrected
  return false;
}

