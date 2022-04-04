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

  // Get the true beam particle
  AnaTrueParticlePD* trueBeamPart = pdAnaUtils::GetTrueBeamParticle(event);
  if (!trueBeamPart) return eventWeight;

  // Get the beam particle identification efficiency. It depends on pdg and nominal beam mom (to be added);
  BinnedParamsParams params;
  int index;
  if(!GetBinValues(abs(trueBeamPart->PDG), trueBeamPart->Momentum, params, index))return eventWeight;  

  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,box.SuccessfulBranch);
   
  // Loop over all TrueParticles in the TPC
  for (Int_t itrue=0;itrue< SystBox->nRelevantTrueObjects; itrue++){      
    //    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(SystBox->RelevantTrueObjects[itrue]);            

    // True-reco association done only once per event in EventWeightBase;
    bool found = (SystBox->RelevantTrueObjectsReco[itrue]!=NULL);

    // compute the weight
    eventWeight *= systUtils::ComputeEffLikeWeight(found, toy, GetIndex(), index, params);
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
