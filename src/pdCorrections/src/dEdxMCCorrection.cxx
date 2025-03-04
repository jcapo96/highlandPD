#include "dEdxMCCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
dEdxMCCorrection::dEdxMCCorrection(){
//********************************************************************

  _params_1 = new BinnedParams(std::string(getenv("PDCORRECTIONSROOT"))+"/data","dEdxMC_1", BinnedParams::k2D_SYMMETRIC,"",true);
  _params_2 = new BinnedParams(std::string(getenv("PDCORRECTIONSROOT"))+"/data","dEdxMC_2", BinnedParams::k2D_SYMMETRIC,"",true);
  _paramsparams_1 = new BinnedParamsParams();
  _paramsparams_2 = new BinnedParamsParams();
  
  _random = new TRandom3();
  _random->SetSeed(1);
}

//********************************************************************
dEdxMCCorrection::~dEdxMCCorrection(){
//********************************************************************

  delete _params_1;
  delete _params_2;
  delete _paramsparams_1;
  delete _paramsparams_2;
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

    //get theta XZ angle
    double thetaXZ = atan(abs(part->DirectionStart[0]/part->DirectionStart[2]))*180/TMath::Pi();
    
    // loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      //the uncorrected dEdx
      double dEdx = part->Hits[2][ihit].dEdx;
      if(dEdx<0)continue;

      //get correction
      Float_t correction   = 1;
      Float_t sigma_landau = 0;
      Float_t sigma_gauss  = 0;
      _params_1->GetInterBinValues(thetaXZ, dEdx, *_paramsparams_1);
      _params_2->GetInterBinValues(thetaXZ, dEdx, *_paramsparams_2);
      correction   = _paramsparams_1->mean;
      sigma_landau = _paramsparams_1->sigma;
      sigma_gauss  = _paramsparams_2->sigma;
      //landau random number generation is not generated over mpv and sigma
      //but over a mean value that is not the same as the mpv.
      //doing numerical experiments I found the relation below to get a 
      //zero mpv for a given sigma.
      double mean_landau  = 0.0067*(1-0.2*sigma_landau)*sigma_landau;
      
      part->Hits[2][ihit].dEdx = dEdx*correction
	+ _random->Landau(mean_landau,sigma_landau)
	+ _random->Gaus(0,sigma_gauss);
    }
  }
}

