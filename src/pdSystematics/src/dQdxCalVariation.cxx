#include "dQdxCalVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
dQdxCalVariation::dQdxCalVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxCal", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  // Initialize Calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize();
 
}

//********************************************************************
void dQdxCalVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

   //Get the systematic source values
  Float_t width;
  GetSigmaValueForBin(0, width); //only 1 bin 

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    if(part->Hits[2].empty())continue;

    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      //get the calibration
      part->Hits[2][ihit].dQdx = part->Hits[2][ihit].dQdx * (1 + width*toy.GetToyVariations(_index)->Variations[0]);
      //recompute dEdx
      _cal->ApplyRecombination(part->Hits[2][ihit]);
    }
    
    //recompute derived quantities
    std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    part->Chi2Proton  = chi2pidprot.first;
    std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    part->Chi2Muon    = chi2pidmuon.first;
    part->Chi2ndf     = chi2pidmuon.second;
  }
}

//********************************************************************
bool dQdxCalVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].dQdx = original->Hits[2][ihit].dQdx;
      part->Hits[2][ihit].dEdx = original->Hits[2][ihit].dEdx;
    }

    // reset the top level quantities
    part->Chi2Proton = original->Chi2Proton;
    part->Chi2ndf    = original->Chi2ndf;
    part->Chi2Muon   = original->Chi2Muon;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dQdxCalVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
