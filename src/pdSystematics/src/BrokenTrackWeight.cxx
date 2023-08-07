#include "BrokenTrackWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "ToyBoxKaon.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
BrokenTrackWeight::BrokenTrackWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","brokenTrack", BinnedParams::k1D_EFF_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h BrokenTrackWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& boxB){
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
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,boxB.SuccessfulBranch); // not really needed
   
  //get bestcandidate
  AnaParticlePD* part = box.Candidates[box.BestCandidateIndex];

  //get parent info to check angle condition
  AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(*static_cast<const AnaEventB*>(&event),part->ParentID));
  if(!parent)return eventWeight;

  //compute cos
  double cos = 0;
  for(int i = 0; i < 3; i++)
    cos += part->DirectionStart[i]*parent->DirectionEnd[i];

  //check brokeness
  bool broken = false;
  if(((part->PositionStart[2]>220 && part->PositionStart[2]<234) ||
      (part->PositionStart[2]>458 && part->PositionStart[2]<466))
     && abs(cos)>0.9){ //broken, contributes to efficiency
    broken = true;
    eventWeight *= systUtils::ComputeEffLikeWeight(broken, toy, GetIndex(), index, params);
  }
  else if((part->PositionStart[2]<220 && part->PositionEnd[2]>220) ||
	  (part->PositionStart[2]<458 && part->PositionEnd[2]>466)){ //not broken, contributes to inefficiency
    broken = false;
    eventWeight *= systUtils::ComputeEffLikeWeight(broken, toy, GetIndex(), index, params);
  }


  return eventWeight;
}

//********************************************************************
bool BrokenTrackWeight::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& recObj) const{
//********************************************************************

  (void)event;      
  return true;
}  
