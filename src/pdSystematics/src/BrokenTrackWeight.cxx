#include "BrokenTrackWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "ToyBoxKaon.hxx"
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

  // Casts the toyBox. The systematic box is not needed
  const ToyBoxKaon& box = *static_cast<const ToyBoxKaon*>(&boxB);   
    
  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,boxB.SuccessfulBranch); // not really needed
   
  //get bestcandidate
  AnaParticlePD* part = box.Candidates[box.BestCandidateIndex];
  //AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  //if(!truePart)return eventWeight;

  //is it in the broken bin?
  if(abs(part->PositionStart[2]-230)<10){
    eventWeight = 0.6;
  }
  else eventWeight = 1.027;
  return eventWeight;
}

//********************************************************************
bool BrokenTrackWeight::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& recObj) const{
//********************************************************************

  (void)event;      
  return true;
}  
