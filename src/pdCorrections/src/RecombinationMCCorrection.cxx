#include "RecombinationMCCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
RecombinationMCCorrection::RecombinationMCCorrection(){
//********************************************************************

  _params = new BinnedParams(std::string(getenv("PDCORRECTIONSROOT"))+"/data","MC_Recombination", BinnedParams::k1D_SYMMETRIC,"",true);

  _cal = new Calorimetry();
  _cal->Initialize();

  Float_t AlphaCorrection = 0;
  Float_t BetaCorrection  = 0;
  AlphaCorrection = _params->GetMeanValueForBin(0, AlphaCorrection); 
  BetaCorrection  = _params->GetMeanValueForBin(1, BetaCorrection); 

  _cal->SetModBoxA(_cal->GetModBoxA()*AlphaCorrection);
  _cal->SetModBoxB(_cal->GetModBoxB()*BetaCorrection);
}

//********************************************************************
RecombinationMCCorrection::~RecombinationMCCorrection(){
//********************************************************************

  delete _params;
  delete _cal;
}

//********************************************************************
void RecombinationMCCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  //this correction is only applied to MC
  if(!spill.GetIsMC())return;

  // Loop over particles
  for(Int_t ipart = 0; ipart < (Int_t)bunch->Particles.size(); ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    //Apply recombination calculation
    _cal->ApplyRecombination(part);
    
    // //recompute derived quantities
    // std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    // part->Chi2Proton  = chi2pidprot.first;
    // std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    // part->Chi2Muon    = chi2pidmuon.first;
    // part->Chi2ndf     = chi2pidmuon.second;
  }
}

