#include "LifetimeVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
LifetimeVariation::LifetimeVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dQdxNorm", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  // Initialize Calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize();
}

//********************************************************************
void LifetimeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //Get the systematic source values
  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin 
 
  //Float_t mc_lt   = _cal->GetLifetime();
  // Float_t data_rQ = 0.89;
  // Float_t data_lt = -2.3/log(data_rQ)*1000; //nominal value used for correction
  // Float_t fake_lt = -2.3/log(data_rQ*(1+sigma*toy.GetToyVariations(_index)->Variations[0]))*1000; //varied value

  Float_t data_lt = 16.54*1000;//-2.3/log(data_rQ)*1000; //nominal value used for correction
  Float_t fake_lt = data_lt*(1+sigma*toy.GetToyVariations(_index)->Variations[0])*1000; //varied value

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    
    if(part->Hits[2].empty())continue;

    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      _cal->SetLifetime(fake_lt);
      _cal->UndoLifetimeCorrection(part->Hits[2][ihit]);
      _cal->SetLifetime(data_lt);
      _cal->ApplyLifetimeCorrection(part->Hits[2][ihit]);
    }
    // //recompute derived quantities
    // std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    // part->Chi2Proton  = chi2pidprot.first;
    // part->Chi2ndf     = chi2pidprot.second;
    // std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    // part->Chi2Muon    = chi2pidmuon.first;
  }
}

//********************************************************************
bool LifetimeVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].dQdx_NoSCE = original->Hits[2][ihit].dQdx_NoSCE;
      part->Hits[2][ihit].dQdx_elife = original->Hits[2][ihit].dQdx_elife;
      //part->Hits[2][ihit].dEdx = original->Hits[2][ihit].dEdx;
    }

    // reset the top level quantities
    // part->Chi2Proton = original->Chi2Proton;
    // part->Chi2ndf    = original->Chi2ndf;
    // part->Chi2Muon   = original->Chi2Muon;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool LifetimeVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
