#include "dQdxCalibVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
dQdxCalibVariation::dQdxCalibVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","dQdxCalib", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void dQdxCalibVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

    if (!original)         continue;
    
    // We need the errors part of the data file but as well the relative uncertainty for sigma
    Float_t mean_corr, mean_var;
    Int_t mean_index;
    
    // Apply the variation to the event model quantities
    // Only plane 2 for the moment
    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){

        AnaHitPD& hit = part->Hits[i][j];                

        // Get the systematics parameters as a function of X
        if (!GetBinValues(hit.Position.x(), mean_corr,  mean_var,  mean_index))  return;
        
        // This is a multiplicative correction assuming the nominal is correct and there is a gaussian error
        // in the multiplicative factors (X_cal, XZ_cal)
        hit.dQdx_NoSCE = hit.dQdx = hit.dQdx*(mean_corr +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);        ;
      }
    }
  }
}

//********************************************************************
bool dQdxCalibVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************


  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=0;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        part->Hits[i][j].dQdx       = original->Hits[i][j].dQdx;
        part->Hits[i][j].dQdx_NoSCE = original->Hits[i][j].dQdx_NoSCE;
      }
    }
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dQdxCalibVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t dQdxCalibVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

