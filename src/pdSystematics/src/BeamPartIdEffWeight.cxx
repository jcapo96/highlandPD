#include "BeamPartIdEffWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
BeamPartIdEffWeight::BeamPartIdEffWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","BeamPartIdEffWeight", BinnedParams::k2D_EFF_ASSYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h BeamPartIdEffWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box){
//********************************************************************
  
  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Get the true beam particle and the event info
  AnaTrueParticlePD* trueBeamPart = pdAnaUtils::GetTrueBeamParticle(event);
  const AnaEventInfoPD* EventInfo = static_cast<const AnaEventInfoPD*>(static_cast<const AnaEventB*>(&event)->EventInfo);
  if(!trueBeamPart || !EventInfo) return eventWeight;

  // Get the beam particle identification efficiency. It depends on pdg and nominal beam mom 
  BinnedParamsParams params;
  int index;
  
  if(!GetBinValues(abs(trueBeamPart->PDG), EventInfo->NominalBeamMom, params, index))return eventWeight;  

  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,box.SuccessfulBranch);

  bool found = false;
  // Loop over all RecoParticles in the TPC
  for(Int_t i = 0; i < SystBox->nRelevantRecObjects; i++){      
    AnaParticlePD* part = static_cast<AnaParticlePD*>(SystBox->RelevantRecObjects[i]);            

    //look for the beam particle
    if(part->isPandora){
      AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
      if(!truePart)break;
      if(truePart->Origin==4)found = true;
      else                   found = false;
      
      // compute the weight
      eventWeight *= systUtils::ComputeEffLikeWeight(found, toy, GetIndex(), index, params);
      break;
    }
  }

  return eventWeight;
}

//********************************************************************
bool BeamPartIdEffWeight::IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const{
//********************************************************************

  (void)event;

  /*const AnaTrueParticleB& truePart = *static_cast<const AnaTrueParticleB*>(&trueObj);

  // consider only charged true particles with momenta below 0.4 (as an example)
  if(truePart.Charge==0) return false;
  if(truePart.Momentum<0.4) return false;*/
    
  return true;
}
