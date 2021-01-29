#include "HitVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
HitVariation::HitVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","Hit", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void HitVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
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
    //for (Int_t j = 0 ; j < (Int_t)NMAXHITSPERPLANE ; j++){
    //  part->HitPosition[2].at(j).SetX(original->HitPosition[2].at(j).X() + mean_corr + (mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);
      //get the corrected coordinates
    //  std::vector<double> point  = {part->HitX[2][j],part->HitY[2][j],part->HitZ[2][j]};
    //  std::vector<double> offset = sce->GetPosOffsets(point);
    //  part->HitX_corrected[2][j] = part->HitX[2][j] - offset[0];
    //}
    
    // Recompute the derived quantities. TODO: This should be done somewhere else, in a more transparent way
    //part->corrected_Length = pdAnaUtils::ComputeCorrectedTrackLength(part,NMAXHITSPERPLANE);
  }
}

//********************************************************************
bool HitVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    //for(Int_t j = 0; j < (Int_t)NMAXHITSPERPLANE; j++){
    //    part->HitX[2][j]           = original->HitX[2][j];
    //    part->HitX_corrected[2][j] = original->HitX_corrected[2][j];
    //     part->corrected_Length     = original->corrected_Length;
    //}   
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool HitVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;

  // True track should always exist
  if (!track.TrueObject) return false;

  AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(track.TrueObject);

  // consider protons, muons and pions
  if(abs(truePart->PDG) == 2212 ) return true;      
  if(abs(truePart->PDG) == 211 )  return true;      
  if(abs(truePart->PDG) == 13 )   return true;      

  return false;
}

//********************************************************************
Int_t HitVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b = 0; b < sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

