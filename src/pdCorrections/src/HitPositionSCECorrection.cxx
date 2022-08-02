#include "HitPositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "Calorimetry.hxx"


//********************************************************************
HitPositionSCECorrection::HitPositionSCECorrection(){
//********************************************************************

  sce = new SpaceCharge();
  sce->Initialize();
  
}

//********************************************************************
void HitPositionSCECorrection::Apply(AnaSpillC& spillC){
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

    for (UInt_t j = 0; j < part->Hits[2].size(); j++){
      TVector3 offset = sce->GetCalPosOffsets(part->Hits[2][j].PositionNoSCE, part->Hits[2][j].TPCid);
      part->Hits[2][j].PositionNoSCE.SetX(part->Hits[2][j].PositionNoSCE.X() - offset.X());
      part->Hits[2][j].PositionNoSCE.SetY(part->Hits[2][j].PositionNoSCE.Y() + offset.Y());
      part->Hits[2][j].PositionNoSCE.SetZ(part->Hits[2][j].PositionNoSCE.Z() + offset.Z());
    }

    //his has nothing to do here but just testing dEdx calibration
    Calorimetry* _cal = new Calorimetry();
    _cal->Initialize();
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      std::cout << "SCE dQdx " << part->Hits[2][ihit].dQdx << " NoSCE dQdx " << part->Hits[2][ihit].dQdx_NoSCE << std::endl;
      _cal->CalibrateHit(part->Hits[2][ihit]);
      std::cout << original->Hits[2][ihit].dEdx << " " << part->Hits[2][ihit].dEdx << std::endl;
    }
  }
}

