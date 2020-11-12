#include "HitPositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
HitPositionSCECorrection::HitPositionSCECorrection(){
//********************************************************************

  
}

//********************************************************************
void HitPositionSCECorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // Loop over all bunches
  for (unsigned int i = 0; i < spill.Bunches.size(); i++){
    AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[i]);

    // Loop over all relevant tracks for this variation
    for (UInt_t itrk = 0; itrk < bunch->Particles.size(); itrk++){
      
      AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[itrk]);
      
      // The un-corrected particle
      const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

      if (!original) continue; //?

      for (UInt_t j = 0; j < part->Hits[2].size(); j++){
        TVector3 offset = sce.GetPosOffsets(part->Hits[2][j].Position);
        part->Hits[2][j].Position = part->Hits[2][j].Position-offset;
      }
    }
  }
}

