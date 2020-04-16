#define pdDataClasses_C

#include "pdDataClasses.hxx"


// define a constant value for uninitialised parameters
//const Float_t  kFloatUnassigned = -999.;
const Double_t kDoubleUnassigned = -999.;
const Int_t    kIntUnassigned = -999;


//********************************************************************
AnaParticlePD::AnaParticlePD():AnaParticle(){
//********************************************************************

  Type = kUnknown;
  isBeamPart    = false;
  isPandoraPart = false;
  
}

//********************************************************************
AnaParticlePD::~AnaParticlePD(){
//********************************************************************

}

//********************************************************************
AnaParticlePD::AnaParticlePD(const AnaParticlePD& part):AnaParticle(part){
//********************************************************************

  Type           = part.Type;
  isBeamPart     = part.isBeamPart;
  isPandoraPart  = part.isPandoraPart;
}

//********************************************************************
void AnaParticlePD::Print() const{
//********************************************************************

  AnaParticle::Print();
  
  std::cout << "-------- AnaParticlePD --------- " << std::endl;
  std::cout << "Type:                    " << Type << std::endl;
  std::cout << "isPandoraPart:           " << isPandoraPart << std::endl;
  std::cout << "PassBeamCut:             " << isBeamPart << std::endl;

  
}


//********************************************************************
AnaTrueParticlePD::AnaTrueParticlePD():AnaTrueParticle(){
//********************************************************************

  Pi0_decay_ID.clear();
  Origin  = kIntUnassigned;
  Matched = false;
  
}

//********************************************************************
AnaTrueParticlePD::~AnaTrueParticlePD(){
//********************************************************************
  
}

//********************************************************************
AnaTrueParticlePD::AnaTrueParticlePD(const AnaTrueParticlePD& truePart):AnaTrueParticle(truePart){
//********************************************************************

  for (UInt_t i=0;i<truePart.Pi0_decay_ID.size();i++)
    Pi0_decay_ID.push_back(truePart.Pi0_decay_ID[i]);
 
  Origin = truePart.Origin;
  Matched = truePart.Matched;

}

//********************************************************************
void AnaTrueParticlePD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrueParticlePD --------- " << std::endl;

  AnaTrueParticle::Print();
  
  if (Pi0_decay_ID.size())
    std::cout << "First Pi0 decay ID:    " << Pi0_decay_ID[0] << std::endl;
  std::cout << "Origin:                " << Origin << std::endl;
  std::cout << "Matched:               " << Matched << std::endl; 
}


//********************************************************************
AnaBeamPD::AnaBeamPD(){
//********************************************************************

  POT           = kIntUnassigned;
  Spill         = kIntUnassigned;

  BeamTrigger   = kIntUnassigned;
  TOF           = kDoubleUnassigned;
  BeamTrackTime = kDoubleUnassigned;
  BeamMomentum  = kDoubleUnassigned;
  BeamMomentumInTPC  = kDoubleUnassigned;
  for (int i=0;i<2;i++){
    CerenkovStatus[i]   = kIntUnassigned;
    CerenkovTime[i]     = kDoubleUnassigned;
    CerenkovPressure[i] = kDoubleUnassigned;
  }

  nMomenta=kIntUnassigned;
  nTracks =kIntUnassigned;
  
  for (int i=0;i<3;i++){
    nFibers[i]=kIntUnassigned;
  }

  PDGs.clear();

}

//********************************************************************
AnaBeamPD::~AnaBeamPD(){
//********************************************************************

}

//********************************************************************
AnaBeamPD::AnaBeamPD(const AnaBeamPD& beam):AnaBeam(beam){
//********************************************************************

  BeamTrigger   = beam.BeamTrigger;
  TOF           = beam.TOF;
  BeamTrackTime = beam.BeamTrackTime;
  BeamMomentum  = beam.BeamMomentum;
  BeamMomentumInTPC  = beam.BeamMomentumInTPC;
  for (int i=0;i<2;i++){
    CerenkovStatus[i]   = beam.CerenkovStatus[i];
    CerenkovTime[i]     = beam.CerenkovTime[i];
    CerenkovPressure[i] = beam.CerenkovPressure[i];
  }

  nMomenta = beam.nMomenta;
  nTracks  = beam.nTracks;
  
  for (int i=0;i<3;i++){
    nFibers[i]=beam.nFibers[i];
  }

  PDGs.clear();
  for (UInt_t i=0; i<beam.PDGs.size(); i++)
    PDGs.push_back(beam.PDGs[i]);
  
}

//********************************************************************
void AnaBeamPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaBeamPD --------- " << std::endl;

  std::cout << "BeamTrigger:      " << BeamTrigger << std::endl;
  std::cout << "TOF:              " << TOF << std::endl;
  std::cout << "BeamTrackTime:    " << BeamTrackTime << std::endl;
  std::cout << "CerenkovStatus:   " << CerenkovStatus[0]   << " " << CerenkovStatus[1]   << std::endl;
  std::cout << "CerenkovTime:     " << CerenkovTime[0]     << " " << CerenkovTime[1]     << std::endl;
  std::cout << "CerenkovPressure: " << CerenkovPressure[0] << " " << CerenkovPressure[1] << std::endl;
  std::cout << "BeamMomentum:     " << BeamMomentum << std::endl;
  std::cout << "BeamMomentumInTPC:" << BeamMomentumInTPC << std::endl;
  std::cout << "nFibers(P1,P2,P3):" << nFibers[0] << " " << nFibers[1] << " " << nFibers[2] << std::endl;
  std::cout << "nMomenta:         " << nMomenta << std::endl;
  std::cout << "nTracks:          " << nTracks << std::endl;

  std::cout << "PDGs:             ";
  for (UInt_t i=0;i<PDGs.size();i++) std::cout << PDGs[i] << " ";
  std::cout << std::endl;

}
