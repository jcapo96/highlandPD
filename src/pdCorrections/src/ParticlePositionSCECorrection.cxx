#include "ParticlePositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
ParticlePositionSCECorrection::ParticlePositionSCECorrection(){
//********************************************************************

  _sce = new SpaceCharge();
  _sce->Initialize();
}

//********************************************************************
void ParticlePositionSCECorrection::Apply(AnaSpillC& spillC){
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

    int ntps = part->TrjPoints.size();
    if(ntps<1)continue;

    //correct trajectory points
    _sce->ApplyTrjPointPositionCorrection(part);
    _sce->ApplyTrjPointDirectionCorrection(part);
    
    //correct length
    part->Length = pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(part);

    //correct end and start position/direction
    pdAnaUtils::ComputeParticlePositionAndDirection(part);
  }
}

