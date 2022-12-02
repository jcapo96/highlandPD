#include "HitPitchSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
HitPitchSCECorrection::HitPitchSCECorrection(){
//********************************************************************

  //initialize space charge effect 
  _sce = new SpaceCharge();
  _sce->Initialize();
  
}

//********************************************************************
void HitPitchSCECorrection::Apply(AnaSpillC& spillC){
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
      double charge = part->Hits[2][j].dQdx_NoSCE * part->Hits[2][j].Pitch_NoSCE;
      TVector3 hitpos = part->Hits[2][j].Position_NoSCE;
      TVector3 hitdir = part->Hits[2][j].Direction_NoSCE;
      
      //compute projection of YZ plane to wire width (Z direction for collection)
      double AngleToVert = 0; //only for collection plane, which is vertical. 
      double cosgamma = abs(sin(AngleToVert)*hitdir.Y() + cos(AngleToVert)*hitdir.Z());
      double pitch = 0.479 / cosgamma; //same values as larsoft with exception of numerical precission
      //std::cout << pitch << " " << part->Hits[2][j].Pitch_NoSCE << std::endl;

      //correct pitch by SCE effect
      TVector3 dirProjection(hitpos.X()+pitch*hitdir.X(),hitpos.Y()+pitch*hitdir.Y(),hitpos.Z()+pitch*hitdir.Z());
      TVector3 dirOffset = _sce->GetCalPosOffsets(dirProjection, part->Hits[2][j].TPCid);
      TVector3 posOffset = _sce->GetCalPosOffsets(hitpos, part->Hits[2][j].TPCid);
      TVector3 dirCorrection(pitch*hitdir.X() - dirOffset.X() + posOffset.X(),
			     pitch*hitdir.Y() + dirOffset.Y() - posOffset.Y(),
			     pitch*hitdir.Z() + dirOffset.Z() - posOffset.Z());
      pitch = dirCorrection.Mag();
      //std::cout << pitch << " " << part->Hits[2][j].Pitch << std::endl;
      std::cout << charge / pitch << " " << part->Hits[2][j].dQdx << std::endl;
    }
  }
}

