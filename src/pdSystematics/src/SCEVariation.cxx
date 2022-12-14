#include "SCEVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCEVariation::SCEVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  //Initialize SCE object to apply variations
  _sce = new SpaceCharge();
  _sce->Initialize();

  // Initialize Calorimetry object using the SCE just created
  _cal = new Calorimetry();
  _cal->SetSCE(_sce);
  _cal->Initialize();
}

//********************************************************************
void SCEVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //Get the systematic source values
  Float_t width;
  GetSigmaValueForBin(0, width); //only 1 bin for the moment

  //Vary SCE map
  _sce->ApplyDisplacementVariation(1 + width*toy.GetToyVariations(_index)->Variations[0]);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);

    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    if(part->Hits[2].empty())continue;
    
    //important, two things are varied: XYZ position of the hit and pitch
    _sce->ApplyPositionCorrection(part);
    _cal->ApplySCECorrection(part);
    _cal->ApplyLifetimeCorrection(part);
  }
}

//********************************************************************
bool SCEVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  _sce->ResetToNominal();
  
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].dQdx = original->Hits[2][ihit].dQdx;
    }
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool SCEVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
