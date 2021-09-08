#include "BeamMomSmearingCorrection.hxx"
#include "pdDataClasses.hxx"
#include <cassert>
#include "Parameters.hxx"

//#define DEBUG

//********************************************************************
BeamMomSmearingCorrection::BeamMomSmearingCorrection(){
//********************************************************************

  _seedValue = ND::params().GetParameterI("protoDuneExampleAnalysis.Corrections.SeedValue");
  // Random generator
  _random.SetSeed(_seedValue);  
}

//********************************************************************
void BeamMomSmearingCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // No correction for data
  if (!spill.GetIsMC())
    return;

  // Get beam
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);

  // Get beam particle
  AnaParticle* beamPart = static_cast<AnaParticle*>(beam->BeamParticle);
  if (!beamPart)                        return; //?
    
  float smearing = 0.04;
      
  beamPart->Momentum = static_cast<AnaTrueParticle*>(beamPart->TrueObject)->Momentum + _random.Gaus(0,smearing);
  beam->BeamMomentum = beamPart->Momentum;
  
  //std::cout << "old momentum " << static_cast<AnaTrueParticle*>(part->TrueObject)->Momentum << " || new momentum " << beam->BeamMomentum << std::endl;

}


