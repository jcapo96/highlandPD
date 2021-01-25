#include "LifetimeVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
LifetimeVariation::LifetimeVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","Lifetime", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void LifetimeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  double fEventT0=500;
  
  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

    if (!original)         continue;
    
    // We need the errors part of the data file but as well the relative uncertainty for sigma
    Float_t mean_corr, mean_var;
    Int_t mean_index;
    
    // Get the systematics parameters
    if (!GetBinValues(0.5, mean_corr,  mean_var,  mean_index))  return;
       
    // Apply the variation to the event model quantities
    // Only plane 2 for the moment
    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){

        AnaHitPD& hit = part->Hits[i][j];                

        // Set the nominal electron lifetime and get the lifetime correction
        fCalo.SetElectronLifetime(mean_corr);
        Float_t corr_nominal = fCalo.LifetimeCorrection(hit.PeakTime, fEventT0);

        // Set the varied lifetime and get the correction
        fCalo.SetElectronLifetime(mean_corr +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);        
        Float_t corr = fCalo.LifetimeCorrection(hit.PeakTime, fEventT0);

        // Uncorrect for the nominal and apply the varied correction only to dQdx
        hit.dQdx_NoSCE = hit.dQdx = hit.dQdx*corr/corr_nominal;
      }
    }
  }
}

//********************************************************************
bool LifetimeVariation::UndoSystematic(AnaEventC& event){
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
bool LifetimeVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t LifetimeVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

