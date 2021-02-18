#include "CNNVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include "HitVariationBase.hxx"

/*
  This is the dEdx calibration variation, applying different variations at the hit level 
  (this avoid redoing the it loop in several places). Those has to be added as daughters 
  of this variation in the analysis algorithm. Some options are:
  - Lifetime
  - dQdxCalib
  - Recombination

  This variation also recomputes all derived quantities depending on dQdx:
  - dEdx
  - Chi2Proton
  - Truncated dEdx
*/

//********************************************************************
CNNVariation::CNNVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","CNN", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void CNNVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // We need the errors part of the data file but as well the relative uncertainty for sigma
  Float_t C_mean, C_var;
  Int_t C_index;

  // Get the systematics parameters
  if (!GetBinValues(0.5,C_mean, C_var, C_index))  return;

  Float_t C = C_mean +  C_var*toy.GetToyVariations(_index)->Variations[C_index]; 

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)         continue;

    // Apply the variation to the event model quantities
    // Only plane 2 for the moment
    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        AnaHitPD& hit = part->Hits[i][j];                
        for (UInt_t k=0;k<hit.Signal.size();k++){
          hit.Signal[k] *= C;
        }
      }
    }

  }

}

//********************************************************************
bool CNNVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        for (UInt_t k=0;k<part->Hits[i][j].Signal.size();k++){
          part->Hits[i][j].Signal[k]       = original->Hits[i][j].Signal[k];
        }
      }
    }
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool CNNVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
  //**************************************************

  (void)event;
  return true;
}
