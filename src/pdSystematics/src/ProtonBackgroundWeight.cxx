#include "ProtonBackgroundWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "ToyBoxKaon.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
ProtonBackgroundWeight::ProtonBackgroundWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","ProtonBackgroundWeight", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h ProtonBackgroundWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& boxB){
//********************************************************************

  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Get the breaking probability for bin 0 (there is a single bin);
  BinnedParamsParams params;
  int index;
  if(!GetBinValues(0.5, params, index))	return eventWeight;

  // Casts the toyBox. The systematic box is not needed
  const ToyBoxKaon& box = *static_cast<const ToyBoxKaon*>(&boxB);   
  if(box.Candidates.empty())return eventWeight;
    
  // Get the SystBox for this event, and the appropriate selection and branch
  //SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,boxB.SuccessfulBranch); // not really needed
   
  //get bestcandidate
  AnaParticlePD* part = box.Candidates[box.BestCandidateIndex];
  if(!part)return eventWeight;

  //get true particle
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return eventWeight;
  
  if(abs(truePart->PDG)!=2212)return eventWeight;
  
  eventWeight.Systematic = params.mean +  params.sigma*toy.GetToyVariations(_index)->Variations[index];
  eventWeight.Correction = params.mean;

  return eventWeight;
}

//********************************************************************
bool ProtonBackgroundWeight::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& recObj) const{
//********************************************************************

  (void)event;      
  return true;
}  
