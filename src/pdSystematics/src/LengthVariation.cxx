#include "LengthVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
LengthVariation::LengthVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","length", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void LengthVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
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
    part->Length = original->Length + (mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);

    // Recompute the derived quantities. TODO: This should be done somewhere else, in a more transparent way
    part->PositionEnd[0] = original->PositionEnd[0]+(part->Length-original->Length)*original->DirectionEnd[0];
    part->PositionEnd[1] = original->PositionEnd[1]+(part->Length-original->Length)*original->DirectionEnd[1];
    part->PositionEnd[2] = original->PositionEnd[2]+(part->Length-original->Length)*original->DirectionEnd[2];
  }
}

//********************************************************************
bool LengthVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    part->Length         = original->Length;
    part->PositionEnd[0] = original->PositionEnd[0];
    part->PositionEnd[1] = original->PositionEnd[1];
    part->PositionEnd[2] = original->PositionEnd[2]; 
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool LengthVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;

  // True track should always exist
  if (!track.TrueObject) return false;

  AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(track.TrueObject);

  // only consider true pions and muons
  if      (abs(truePart->PDG) == 211 ) return true;      
  else if (abs(truePart->PDG) == 13)   return true;

  return false;
}

//********************************************************************
Int_t LengthVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

