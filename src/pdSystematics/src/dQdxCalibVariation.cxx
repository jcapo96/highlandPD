#include "dQdxCalibVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
dQdxCalibVariation::dQdxCalibVariation(): HitVariationBase(){
//********************************************************************

  _paramX    = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCalibX",  BinnedParams::k1D_SYMMETRIC);
  _paramYZ   = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCalibYZ", BinnedParams::k1D_SYMMETRIC);
  _paramNQ   = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxNQ",      BinnedParams::k1D_SYMMETRIC);
  _paramCcal = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCcal",    BinnedParams::k1D_SYMMETRIC);
  
  // Read the systematic source parameters from the data files
  SetNParameters(_paramX->GetNBins()+_paramYZ->GetNBins()+_paramNQ->GetNBins()+_paramCcal->GetNBins());
}

//********************************************************************
void dQdxCalibVariation::Apply(const ToyExperiment& toy, AnaHitPD& hit){
//********************************************************************

  //--------- If there is a single bin this can be done in constructor ----
  // Get the systematics parameters
  Float_t X_mean, X_var, YZ_mean, YZ_var, NQ_mean, NQ_var, Ccal_mean, Ccal_var;
  Int_t X_index, YZ_index, NQ_index, Ccal_index;

  if (   !_paramX->GetBinValues(0.5, X_mean,    X_var,    X_index))    return;
  if (  !_paramYZ->GetBinValues(0.5, YZ_mean,   YZ_var,   YZ_index))   return;
  if (  !_paramNQ->GetBinValues(0.5, NQ_mean,   NQ_var,   NQ_index))   return;
  if (!_paramCcal->GetBinValues(0.5, Ccal_mean, Ccal_var, Ccal_index)) return;
  //-----------------------------------------------------------------------

  //--------- If there is a single bin this can be done in Apply(toy, event) ----
  // do the variations
  Float_t C_X  =  X_mean   +     X_var*toy.GetToyVariations(_index)->Variations[X_index];
  Float_t C_YZ = YZ_mean   +    YZ_var*toy.GetToyVariations(_index)->Variations[X_index+_paramX->GetNBins()];
  Float_t NQ   = NQ_mean   +    NQ_var*toy.GetToyVariations(_index)->Variations[X_index+_paramX->GetNBins()+_paramYZ->GetNBins()];
  Float_t Ccal = Ccal_mean +  Ccal_var*toy.GetToyVariations(_index)->Variations[X_index+_paramX->GetNBins()+_paramYZ->GetNBins()+_paramNQ->GetNBins()];
  //-----------------------------------------------------------------------------
  
  hit.dQdx       *= C_X*C_YZ*NQ*Ccal;
  hit.dQdx_NoSCE *= C_X*C_YZ*NQ*Ccal;

}

//********************************************************************
bool dQdxCalibVariation::UndoSystematic(const AnaHitPD& original_hit, AnaHitPD& hit){
  //********************************************************************
  
  hit.dQdx       = original_hit.dQdx;
  hit.dQdx_NoSCE = original_hit.dQdx_NoSCE;

  // Don't reset the spill to corrected
  return false;
}

