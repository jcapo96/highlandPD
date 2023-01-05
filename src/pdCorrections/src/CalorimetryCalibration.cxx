#include "CalorimetryCalibration.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
CalorimetryCalibration::CalorimetryCalibration(){
//********************************************************************

  //initialize space charge effect
  _sce = new SpaceCharge();
  _sce->Initialize();
  //initialize calorimetry object
  _cal = new Calorimetry();
  _cal->Initialize(_sce); //use same sce to save memory
}

//********************************************************************
void CalorimetryCalibration::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  // Loop over particles
  for(UInt_t ipart = 0; ipart < bunch->Particles.size(); ipart++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    //loop over hits
    for(UInt_t ihit = 0; ihit < part->Hits[2].size(); ihit++){
      _cal->ApplySCECorrection(part->Hits[2][ihit]);      //pitch by SCE
      _cal->ApplyLifetimeCorrection(part->Hits[2][ihit]); //dQdx by lifetime
      _cal->CalibratedQdx(part->Hits[2][ihit]);           //dQdx XYZ
      _cal->ApplyRecombination(part->Hits[2][ihit]);      //dEdx by recombination
      //it would be nice to think appropiate names and doing it in a 
      //centralized and more comfortable way
    }
  }
}

