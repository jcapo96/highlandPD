#include "LifetimeVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
LifetimeVariation::LifetimeVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","Lifetime", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  // Initialize Calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize();
}
//********************************************************************
LifetimeVariation::~LifetimeVariation(){
//********************************************************************

  if(_cal)delete _cal;
}

//********************************************************************
void LifetimeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);
  const AnaEventPD* eventpd = static_cast<const AnaEventPD*>(&event);

  //Get the systematic source values
  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin 
  
  //set appropiate lifetime depending on mc/data and run
  _cal->SetLifetime(*eventpd);

  //compute varied lt from varied rQ
  Float_t nom_lt = _cal->GetLifetime()/1000; //us to ms
  Float_t nom_rQ = exp(-2.3/nom_lt);
  Float_t var_lt = -2.3/log(nom_rQ*(1+sigma*toy.GetToyVariations(_index)->Variations[0])); //varied lt

  //set new lifetime
  _cal->SetLifetime(var_lt*1000); //ms to us

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    
    if(part->Hits[2].empty())continue;

    //correct every hit 
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
      _cal->ApplyLifetimeCorrection(part->Hits[2][ihit]);
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
      part->Hits[2][ihit].dQdx_SCE   = original->Hits[2][ihit].dQdx_SCE;
      part->Hits[2][ihit].dQdx_elife = original->Hits[2][ihit].dQdx_elife;
    }
  }

  //reset calorimetry object
  _cal->ResetLifetime();
  
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
