#include "BeamMomCorrection.hxx"
#include "pdDataClasses.hxx"
#include <cassert>

//#define DEBUG

//********************************************************************
BeamMomCorrection::BeamMomCorrection(){
//********************************************************************

}

//********************************************************************
void BeamMomCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // No correction for MC
  if (spill.GetIsMC())
    return;

  // Get beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticle* beamPart = static_cast<AnaParticle*>(beam->BeamParticle);
  if (!beamPart) return;
    
  // The un-corrected particle
  const AnaParticle* original = static_cast<const AnaParticle*>(beamPart->Original); 
  if (!original) return; 
  
  float corr = 1.041883726;
      
  beamPart->Momentum = original->Momentum * corr;
  beam->BeamMomentum = beam->BeamMomentum * corr;
  
  //  std::cout << "old momentum " << original->Momentum << " || new momentum " << beamPart->Momentum << std::endl;
}


