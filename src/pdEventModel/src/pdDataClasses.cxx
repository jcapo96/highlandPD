#define pdDataClasses_C

#include "pdDataClasses.hxx"


// define a constant value for uninitialised parameters
//const Float_t  kFloatUnassigned = -999.;
const Double_t kDoubleUnassigned = -999.;
const Int_t    kIntUnassigned = -999;
const Float_t  kFloatUnassigned = -999.;

//********************************************************************
AnaParticlePD::AnaParticlePD():AnaParticle(){
//********************************************************************

  Type = kUnknown;
  isBeamPart    = false;
  isPandora = false;


  FitPDG        = kFloatUnassigned;
    for (Int_t i=0;i<3;i++){
    NHitsPerPlane[i] = kIntUnassigned;
    for (UInt_t j=0;j<NMAXHITSPERPLANE;j++){
      dEdx[i][j]=kFloatUnassigned;
      dQdx[i][j]=kFloatUnassigned;
      dEdx_corr[i][j]=kFloatUnassigned;
      dQdx_corr[i][j]=kFloatUnassigned;
      HitX[i][j]=kFloatUnassigned;
      HitY[i][j]=kFloatUnassigned;
      HitZ[i][j]=kFloatUnassigned;
      ResidualRange[i][j]=kFloatUnassigned;
    }
  }  



  for (int i=0; i<3; i++) {
    PIDA[i]=kFloatUnassigned;
    ReconPDG[i]=kIntUnassigned;
    for (int j=0; j<10; j++) {
      PID[i][j]=kFloatUnassigned;
      CALO[i][j]=kFloatUnassigned;
    }
  }


  Chi2Proton=kFloatUnassigned;
  Chi2ndf=kIntUnassigned;

  for (int i=0; i<3; i++)
    CNNscore[i]=kFloatUnassigned;
  
  for (int i=0; i<2; i++)
    RangeMomentum[i] = kFloatUnassigned;    
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
  isPandora  = part.isPandora;

  FitPDG        = part.FitPDG;

  
  for (Int_t i=0;i<3;i++){
    NHitsPerPlane[i] = part.NHitsPerPlane[i];
    for (Int_t j=0;j<std::min((Int_t)NMAXHITSPERPLANE,NHitsPerPlane[i]);j++){
      dEdx[i][j]=part.dEdx[i][j];
      dQdx[i][j]=part.dQdx[i][j];
      dEdx_corr[i][j]=part.dEdx_corr[i][j];
      dQdx_corr[i][j]=part.dQdx_corr[i][j];
      HitX[i][j]=part.HitX[i][j];
      HitY[i][j]=part.HitY[i][j];
      HitZ[i][j]=part.HitZ[i][j];
      ResidualRange[i][j]=part.ResidualRange[i][j];
    }
  }



  for (int i=0; i<3; i++) {
    PIDA[i]=part.PIDA[i];
    ReconPDG[i]=part.ReconPDG[i];

    for (int j=0; j<10; j++) {
      PID[i][j]=part.PID[i][j];
      CALO[i][j]=part.CALO[i][j];
    }
  }

  Chi2Proton = part.Chi2Proton;
  Chi2ndf    = part.Chi2ndf;

  for (int i=0; i<3; i++)
    CNNscore[i]=part.CNNscore[i];
  
  
  for (int i=0; i<2; i++)
    RangeMomentum[i] = part.RangeMomentum[i];  
}

//********************************************************************
void AnaParticlePD::Print() const{
//********************************************************************

  AnaParticle::Print();
  
  std::cout << "-------- AnaParticlePD --------- " << std::endl;
  std::cout << "Type:                    " << Type << std::endl;
  std::cout << "isPandora:           " << isPandora << std::endl;
  std::cout << "PassBeamCut:             " << isBeamPart << std::endl;

  std::cout << "FitPDG:                  " << FitPDG << std::endl;
  
  std::cout << "PIDA:                    ";
  for (int i=0;i<3;i++) std::cout << PIDA[i] << " ";
  std::cout << std::endl;

  std::cout << "ReconPDG:                ";
  for (int i=0;i<3;i++) std::cout << ReconPDG[i] << " ";
  std::cout << std::endl;
 
  std::cout << "PID[0]:                  ";
  for (int i=0;i<10;i++) std::cout << PID[0][i] << " ";
  std::cout << std::endl;

  std::cout << "PID[1]:                  ";
  for (int i=0;i<10;i++) std::cout << PID[1][i] << " ";
  std::cout << std::endl;

  std::cout << "PID[2]:                  ";
  for (int i=0;i<10;i++) std::cout << PID[2][i] << " ";
  std::cout << std::endl;

  std::cout << "CALO[0]:                  ";
  for (int i=0;i<10;i++) std::cout << CALO[0][i] << " ";
  std::cout << std::endl;

  std::cout << "CALO[1]:                  ";
  for (int i=0;i<10;i++) std::cout << CALO[1][i] << " ";
  std::cout << std::endl;


  std::cout << "CALO[2]:                  ";
  for (int i=0;i<10;i++) std::cout << CALO[2][i] << " ";
  std::cout << std::endl;

  std::cout << "Chi2Proton:              " << Chi2Proton << std::endl;
  std::cout << "Chi2ndf:                 " << Chi2ndf << std::endl;

  std::cout << "CNN score:               " << CNNscore[0] << " " << CNNscore[1] << " " << CNNscore[2] << std::endl;
  
  std::cout << "RangeMomentum            " << RangeMomentum[0] << " " << RangeMomentum[1] << std::endl;

    
}


//********************************************************************
AnaTrueParticlePD::AnaTrueParticlePD():AnaTrueParticle(){
//********************************************************************

  Pi0_decay_ID.clear();
  Origin  = kIntUnassigned;
  Matched = false;
  LengthInTPC   = kFloatUnassigned;
  MomentumInTPC = kFloatUnassigned;
  
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
  LengthInTPC = truePart.LengthInTPC;
  MomentumInTPC = truePart.MomentumInTPC;
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
  std::cout << "LengthInTPC:           " << LengthInTPC << std::endl;
  std::cout << "MomentumInTPC:         " << MomentumInTPC << std::endl;

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
