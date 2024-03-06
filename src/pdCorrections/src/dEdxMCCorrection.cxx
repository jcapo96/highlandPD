#include "dEdxMCCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
dEdxMCCorrection::dEdxMCCorrection(){
//********************************************************************

  _params = new BinnedParams(std::string(getenv("PDCORRECTIONSROOT"))+"/data","dEdxMC", BinnedParams::k1D_SYMMETRIC,"",true);

}

//********************************************************************
void dEdxMCCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  //if not MC, return
  if(!spill.GetIsMC())return;

  // Loop over particles
  for(UInt_t ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
    
    //get particle 
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    //get correction
    Float_t correction = 1;
    _params->GetMeanValueForBin(0, correction); //only 1 bin
    
    // loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      //the uncorrected dEdx
      double dEdx = part->Hits[2][ihit].dEdx;
      if(dEdx<0)continue;
      part->Hits[2][ihit].dEdx = dEdx*correction;   
    }
  }
}

