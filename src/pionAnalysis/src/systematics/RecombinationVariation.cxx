#include "RecombinationVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"


// MC:   ModBoxA = 0.937±0.020 (ArgoNeuT: 0.93±0.02), ModBoxB = 0.209±0.006 (ArgoNeuT: 0.212±0.002)
// DATA: ModBoxA = 0.854±0.017 (ArgoNeuT: 0.93±0.02) ,ModBoxB = 0.208±0.006 (ArgoNeuT: 0.212±0.002)
                                                         

//********************************************************************
RecombinationVariation::RecombinationVariation(): EventVariationBase(){
//********************************************************************



  _modBoxA = new BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","RecombinationA", BinnedParams::k1D_SYMMETRIC);
  _modBoxB = new BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","RecombinationB", BinnedParams::k1D_SYMMETRIC);
  
  // Read the systematic source parameters from the data files
  SetNParameters(_modBoxA->GetNBins()+_modBoxB->GetNBins());
}

//********************************************************************
void RecombinationVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
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
    Float_t A_mean, A_var, B_mean, B_var;
    Int_t A_index, B_index;

    
    // Get the systematics parameters as a function of X
    if (!_modBoxA->GetBinValues(0.5, A_mean,  A_var,  A_index))  return;
    if (!_modBoxB->GetBinValues(0.5, B_mean,  B_var,  B_index))  return;

    Float_t ModBoxA = A_mean +  A_var*toy.GetToyVariations(_index)->Variations[A_index];
    Float_t ModBoxB = B_mean +  B_var*toy.GetToyVariations(_index)->Variations[B_index]; 
    
    // Set the varied parameters
    fCalo.SetModBoxA(ModBoxA);
    fCalo.SetModBoxB(ModBoxB);        
    
    // Apply the variation to the event model quantities
    // Only plane 2 for the moment
    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){

        AnaHitPD& hit = part->Hits[i][j];                
        
        // compute the dEdx from the varied dQdx
        hit.dEdx      = fCalo.dEdx_from_dQdx(hit.dQdx);
        hit.dEdx_corr = fCalo.dEdx_from_dQdx(hit.dQdx_corr);
      }
    }
  }
}

//********************************************************************
bool RecombinationVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************


  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=0;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        part->Hits[i][j].dEdx      = original->Hits[i][j].dEdx;
        part->Hits[i][j].dEdx_corr = original->Hits[i][j].dEdx_corr;
      }
    }
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool RecombinationVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t RecombinationVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

