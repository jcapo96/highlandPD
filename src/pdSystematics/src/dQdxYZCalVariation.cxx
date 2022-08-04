#include "dQdxYZCalVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
dQdxYZCalVariation::dQdxYZCalVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxYZCal", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  // Initialize Calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize();
}

//********************************************************************
void dQdxYZCalVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    if(part->Hits[2].empty())continue;

    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //Get the systematic source values
    Float_t width;
    GetSigmaValueForBin(0, width); //only 1 bin 
    std::cout << "BIN " << width << std::endl;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].PlaneID = 2; //not filled in the event model, fill it by hand for the moment
      //here it is assumed calibrated dQdx is contained in the event model
      //get the correction
      double YZCal = _cal->GetYZCalibration(part->Hits[2][ihit]);
      //undo the correction and apply the varied correction
      //std::cout << "DQDYZ RATIO " << (YZCal + width*toy.GetToyVariations(_index)->Variations[0]) / YZCal << std::endl;
      part->Hits[2][ihit].dQdx = part->Hits[2][ihit].dQdx * (YZCal + width*toy.GetToyVariations(_index)->Variations[0]) / YZCal;
      //recompute dEdx
      _cal->ApplyRecombination(part->Hits[2][ihit]);
      //std::cout << "DEDYZ RATIO " << part->Hits[2][ihit].dEdx / original->Hits[2][ihit].dEdx << std::endl;
    }
    
    //recompute derived quantities
    std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    part->Chi2Proton  = chi2pidprot.first;
    part->Chi2ndf     = chi2pidprot.second;
    std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    part->Chi2Muon    = chi2pidmuon.first;
  }
}

//********************************************************************
bool dQdxYZCalVariation::UndoSystematic(AnaEventC& event){
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
bool dQdxYZCalVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
