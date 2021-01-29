#include "dEdxVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
dEdxVariation::dEdxVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","dEdx", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void dEdxVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

    if (!part->TrueObject) continue;
    if (!original)         continue;

    AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(part->TrueObject);
    
    // We need the errors part of the data file but as well the relative uncertainty for sigma
    Float_t mean_corr, mean_var;
    Int_t mean_index;
    
    // Get the systematics parameters for this particle type
    if (!GetBinValues(abs(truePart->PDG), mean_corr,  mean_var,  mean_index))  return;

    // Apply the variation to the evnt model quantities
    for (Int_t i=0;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        part->Hits[i][i].dEdx      = original->Hits[i][i].dEdx      *(1 +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]/mean_corr);
        part->Hits[i][i].dEdx_corr = original->Hits[i][i].dEdx_corr *(1 +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]/mean_corr);
        // A variation on the residual range and the track length is also expected. Should this be done here or in a different systematic ?
      }
    }
  }
}

//********************************************************************
bool dEdxVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=0;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        part->Hits[i][i].dEdx      = original->Hits[i][i].dEdx;
        part->Hits[i][i].dEdx_corr = original->Hits[i][i].dEdx_corr;
      }
    }
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dEdxVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;

  // True track should always exist
  if (!track.TrueObject) return false;

  AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(track.TrueObject);

  // only consider true protons, pions, muons and electrons
  if      (abs(truePart->PDG) == 211 ) return true;      
  else if (abs(truePart->PDG) == 2212) return true;
  else if (abs(truePart->PDG) == 13)   return true;
  else if (abs(truePart->PDG) == 11)   return true;
  else if (abs(truePart->PDG) == 321)  return true;

  return false;
}

//********************************************************************
Int_t dEdxVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

