#include "ParticlePositionSCE.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
ParticlePositionSCE::ParticlePositionSCE(){
//********************************************************************

  _sce = new SpaceCharge();
  _sce->Initialize();
}

//********************************************************************
void ParticlePositionSCE::Apply(AnaSpillC& spillC){
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

    TVector3 pos;
    int TPCid = 0;

    //position start
    pos.SetXYX(part->PositionStart[0],part->PositionStart[1],part->PositionStart[2]);
    if(pos.X() < 0)TPCid = 1;
    else           TPCid = 2;
    TVector3 offset = _sce->GetCalPosOffset(pos,TPCid);
    part->PositionStart[0] = pos.X() - offset.X();
    part->PositionStart[1] = pos.Y() + offset.Y();
    part->PositionStart[2] = pos.Z() + offset.Z();

    //position end
    pos.SetXYX(part->PositionEnd[0],part->PositionEnd[1],part->PositionEnd[2]);
    if(pos.X() < 0)TPCid = 1;
    else           TPCid = 2;
    TVector3 offset = _sce->GetCalPosOffset(pos,TPCid);
    part->PositionEnd[0] = pos.X() - offset.X();
    part->PositionEnd[1] = pos.Y() + offset.Y();
    part->PositionEnd[2] = pos.Z() + offset.Z();
  }
}

