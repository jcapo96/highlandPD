#include "BeamMomSmearingCorrection.hxx"
#include "pdDataClasses.hxx"
#include <cassert>
#include "Parameters.hxx"

//********************************************************************
BeamMomSmearingCorrection::BeamMomSmearingCorrection(){
//********************************************************************

  _seedValue = ND::params().GetParameterI("pdCorrections.BeamMomSmearing.SeedValue");
  _NominalBeamMom = ND::params().GetParameterI("pdCorrections.BeamMomSmearing.NominalBeamMom");
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
    
  float smearing = 0;
  if(_NominalBeamMom == 1)smearing = 0.025;
  else if(_NominalBeamMom == 2)smearing = 0.124;
  else if(_NominalBeamMom == 3)smearing = 0.210;
  else{
    std::cout << "beam mom " << _NominalBeamMom << " not contemplated" << std::endl;
    std::exit(1);
  }

  //in PDSPAnalyzer ntuples the beam momentum appears centered in one always
  beamPart->Momentum = beamPart->Momentum*_NominalBeamMom + _random.Gaus(0,smearing);
  beam->BeamMomentum = beamPart->Momentum;
}


