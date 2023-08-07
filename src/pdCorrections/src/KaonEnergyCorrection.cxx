#include "KaonEnergyCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
KaonEnergyCorrection::KaonEnergyCorrection(){
//********************************************************************

  _params = new BinnedParams(std::string(getenv("PDCORRECTIONSROOT"))+"/data","KaonEnergy", BinnedParams::k1D_SYMMETRIC,"",true);

}

//********************************************************************
void KaonEnergyCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  //if not MC, return
  if(!spill.GetIsMC())return;

  BinnedParamsParams par;

  // Loop over particles
  for(UInt_t ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
    
    //get particle and check it is a kaon
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
    if(!truePart)continue;
    if(abs(truePart->PDG)!=321)continue;
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    // loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      //the uncorrected dEdx
      double dEdx = part->Hits[2][ihit].dEdx;
      if(dEdx<0)continue;
      _params->GetInterBinValues(dEdx, par);
      Float_t correction = par.mean;
      part->Hits[2][ihit].dEdx = dEdx*correction;   
    }
  }
}

