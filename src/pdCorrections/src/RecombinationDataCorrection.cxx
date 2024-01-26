#include "RecombinationDataCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
RecombinationDataCorrection::RecombinationDataCorrection(){
//********************************************************************

}

//********************************************************************
void RecombinationDataCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  //this correction is only applied to data
  if(spill.GetIsMC())return;

  //electric field and lar density
  double Efield = 0.553;
  double rho = 1.38877;

  //recombination parameters
  double alpha = 0.930;
  double beta  = 0.212;
  double alpha_new = 0.905;
  double beta_new  = 0.220;

  //calibration constant ratio
  double C_ratio = 0.95;

  // Loop over particles
  for(Int_t ipart = 0; ipart < (Int_t)bunch->Particles.size(); ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    //loop over hits
    for(Int_t ihit = 0; ihit < (Int_t)part->Hits[2].size(); ihit++){
      double dEdx = original->Hits[2][ihit].dEdx;
      part->Hits[2][ihit].dEdx = (C_ratio * exp(beta_new/beta * log(alpha + (beta*dEdx)/(rho*Efield) )) - alpha_new) * (rho*Efield/beta_new);
    }
  }
}

