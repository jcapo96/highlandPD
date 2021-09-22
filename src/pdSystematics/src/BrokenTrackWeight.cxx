#include "BrokenTrackWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "ToyBoxPD.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
BrokenTrackWeight::BrokenTrackWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","brokenTrack", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h BrokenTrackWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& boxB){
//********************************************************************

  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Systematic source parameters
  BinnedParamsParams params;
  int index;

  // Casts the toyBox
  const ToyBoxPD& box = *static_cast<const ToyBoxPD*>(&boxB);   
    
  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,boxB.SuccessfulBranch);
   
  // Loop over all relevant objects
  for (Int_t irec=0;irec< SystBox->nRelevantRecObjects; irec++){      

    AnaParticlePD* part = static_cast<AnaParticlePD*>(SystBox->RelevantRecObjects[irec]);            
    if (!part) continue;

    // use the candidate only for the moment
    if (part!=box.MainTrack) continue;

    // Get the systematic source parameters for this track (depend on theta_y)
    if(!GetBinValues(part->DirectionEnd[1], params, index))	return eventWeight;
    
    // retrieve the end position in Z
    Float_t endpos_z = part->PositionEnd[2];

    // Broken track found
    bool found = (endpos_z < 234);

    // compute the weight
    eventWeight *= systUtils::ComputeEffLikeWeight(found, toy, GetIndex(), index, params);
  }

  return eventWeight;
}

//********************************************************************
bool BrokenTrackWeight::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& recObj) const{
//********************************************************************

  (void)event;
  
  const AnaParticleB& part = *static_cast<const AnaParticleB*>(&recObj);

  // retrieve the end position in Z
  Float_t endpos_z = part.PositionEnd[2];

  // consider only tracks with end z > 222
  if (endpos_z < 222) return false;
      
  return true;
}


/*

  
  (void)event;

  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Gets the selected track from the box
  const ToyBoxPD& box = *static_cast<const ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return eventWeight;
  
  // retrieve the end position in Z
  Float_t endpos_z = box.MainTrack->PositionEnd[2];

  // nothing to do for tracks in first APA upstream of the conflictive region
  if (endpos_z < 222) return eventWeight;

  // fraction of tracks above 222 which end in the 222-234 region
  Float_t broken_prob_data;
  Float_t broken_prob_mc;   
  Float_t c_error;
  Int_t weight_index, weight_index_dummy;

  // Get the systematics parameters 
  if (!GetBinValues(0, broken_prob_data,  c_error,  weight_index))        return eventWeight;
  if (!GetBinValues(1, broken_prob_mc,    c_error,  weight_index_dummy))  return eventWeight;
  
  // nominal data/mc ratio
  Float_t c0 = broken_prob_data/broken_prob_mc; 

  // varied data/mc ratio
  Float_t c  = c0 + c_error*toy.GetToyVariations(_index)->Variations[weight_index];

  // compute the weight
  if (endpos_z < 234){     
    eventWeight.Systematic = c;
    eventWeight.Correction = c0;
  }
  else{
    eventWeight.Systematic = (1 - c *broken_prob_data)/(1 - broken_prob_mc);
    eventWeight.Correction = (1 - c0*broken_prob_data)/(1 - broken_prob_mc);
  }

  return eventWeight;
}
*/  
