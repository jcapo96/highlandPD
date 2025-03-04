#include "RecombinationVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
RecombinationVariation::RecombinationVariation():EventVariationBase(),BinnedParams(){
//********************************************************************

  _BP_A = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","RecombinationA", BinnedParams::k1D_SYMMETRIC);
  _BP_B = new BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","RecombinationB", BinnedParams::k1D_SYMMETRIC);

  // Read the systematic source parameters from the data files
  SetNParameters(_BP_A->GetNBins()+_BP_B->GetNBins());

  // Initialize Calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize();
}

//********************************************************************
void RecombinationVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //Get mean value and relative uncertainty
  Float_t A_mean, A_var, B_mean, B_var;
  Int_t A_index, B_index;
  
  // Get the systematics parameters (This can be outside the loop if it does not depend on Particle properties)
  if(!_BP_A->GetBinValues(0.5, A_mean,  A_var,  A_index))return;
  if(!_BP_B->GetBinValues(0.5, B_mean,  B_var,  B_index))return;

  // They are correlated. To throw correlated random numbers, we have to generate a vector of 2 random numbers and 
  // multiply it by the Cholesky decomposition of the correlation matrix.
  double cholesky[2][2] = {{1,0},{0.85,1}};

  // these are the uncorrelated random numbers
  double var[2] = {toy.GetToyVariations(_index)->Variations[A_index],
		   toy.GetToyVariations(_index)->Variations[B_index+_BP_A->GetNBins()]};

  // now we compute the correlated random numbers
  double cvar[2] = {0};
  for(int irow = 0; irow < 2; irow++){
    double product = 0;
    for(int icolumn = 0; icolumn < 2; icolumn++)
      product += cholesky[irow][icolumn]*var[icolumn];
    cvar[irow] = product;
  }

  // and now we modify the modbox model parameters in a correlated way
  Float_t ModBoxA = _cal->GetModBoxA()*(A_mean +  A_var*cvar[0]);
  Float_t ModBoxB = _cal->GetModBoxB()*(B_mean +  B_var*cvar[1]); 

  //modify alpha and beta values on calorimetry class
  _cal->SetModBoxA(ModBoxA);
  _cal->SetModBoxB(ModBoxB);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);

    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //Apply recombination calculation
    _cal->ApplyRecombination(part);
    
    //recompute derived quantities
    std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    part->Chi2Proton  = chi2pidprot.first;
    std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    part->Chi2Muon    = chi2pidmuon.first;
    part->Chi2ndf     = chi2pidmuon.second;
  }

  _cal->ResetModBoxParameters(); //not really needed but just in case
}

//********************************************************************
bool RecombinationVariation::UndoSystematic(AnaEventC& event){
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
  
  _cal->ResetModBoxParameters(); //double check

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool RecombinationVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
