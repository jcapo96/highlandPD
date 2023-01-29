#define pdDataClasses_C

#include "pdDataClasses.hxx"
#include "AnalysisUtils.hxx"


// define a constant value for uninitialised parameters
//const Float_t  kFloatUnassigned = -999.;
const Double_t kDoubleUnassigned = -999.;
const Int_t    kIntUnassigned = -999;
const Float_t  kFloatUnassigned = -999.;


//********************************************************************
AnaWireID::AnaWireID(const AnaWireID& wireID){
//********************************************************************

  Wire       = wireID.Wire;
  Plane      = wireID.Plane;
  TPC        = wireID.TPC;
  Cryostat  = wireID.Cryostat;
}

//********************************************************************
AnaHitPD::AnaHitPD(){
//********************************************************************

  WireID.Plane  = kIntUnassigned; 
  Integral      = kFloatUnassigned;     
  PeakTime      = kFloatUnassigned;     
  PeakAmplitude = kFloatUnassigned;
  TPCid         = kIntUnassigned;
  PlaneID       = kIntUnassigned;
  Channel       = kIntUnassigned;
  View          = kIntUnassigned; 

  StartTick     = (UInt_t)kIntUnassigned;
  EndTick       = (UInt_t)kIntUnassigned; 

  dQdx_NoSCE          = kFloatUnassigned;
  dEdx_NoSCE          = kFloatUnassigned;
  ResidualRange_NoSCE = kFloatUnassigned;
  Pitch_NoSCE         = kFloatUnassigned;
  Position_NoSCE      = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
  Direction_NoSCE     = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);

  dQdx_SCE = kFloatUnassigned;
  dEdx_SCE = kFloatUnassigned;
  ResidualRange_SCE = kFloatUnassigned;

  dQdx_elife = kFloatUnassigned;
  dEdx_elife = kFloatUnassigned;
  ResidualRange_elife = kFloatUnassigned;

  dQdx          = kFloatUnassigned;
  dEdx          = kFloatUnassigned;
  ResidualRange = kFloatUnassigned;
  Pitch         = kFloatUnassigned;
  Position      = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
  Direction     = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);

  dEdx_calib    = kFloatUnassigned;
  
  Signal.clear();
  //CNN.resize(3);
  CNN[0]=CNN[1]=CNN[2]=kFloatUnassigned;
}

//********************************************************************
AnaHitPD::AnaHitPD(const AnaHitPD& hit){
//********************************************************************

  WireID.Plane  = hit.WireID.Plane;
  WireID.Wire   = hit.WireID.Wire; 
  Integral      = hit.Integral;     
  PeakTime      = hit.PeakTime;     
  PeakAmplitude = hit.PeakAmplitude;
  TPCid         = hit.TPCid;
  PlaneID       = hit.PlaneID;
  Channel       = hit.Channel;
  View          = hit.View;
  Signal        = hit.Signal;
  //CNN.resize(3);
  for (size_t i=0;i<3;i++)
    CNN[i]           = hit.CNN[i];

  StartTick     = hit.StartTick;
  EndTick       = hit.EndTick  ; 
  
  dQdx_NoSCE          = hit.dQdx_NoSCE;
  dEdx_NoSCE          = hit.dEdx_NoSCE;
  ResidualRange_NoSCE = hit.ResidualRange_NoSCE;
  Pitch_NoSCE         = hit.Pitch_NoSCE;
  Position_NoSCE      = hit.Position_NoSCE;
  Direction_NoSCE     = hit.Direction_NoSCE;

  dQdx_SCE = hit.dQdx_SCE;
  dEdx_SCE = hit.dEdx_SCE;
  ResidualRange_SCE = hit.ResidualRange_SCE;

  dQdx_elife = hit.dQdx_elife;
  dEdx_elife = hit.dEdx_elife;
  ResidualRange_elife = hit.ResidualRange_elife;

  dQdx          = hit.dQdx;
  dEdx          = hit.dEdx;
  ResidualRange = hit.ResidualRange;
  Pitch         = hit.Pitch;
  Position      = hit.Position;
  Direction     = hit.Direction;

  dEdx_calib     = hit.dEdx_calib;
}

//********************************************************************
void AnaHitPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaHitPD --------- " << std::endl;

  std::cout << "WireID.Plane:  " << WireID.Plane  << std::endl;
  std::cout << "WireID.Wire:   " << WireID.Wire  << std::endl;
  std::cout << "Channel:       " << Channel  << std::endl;
  std::cout << "View:          " << View << std::endl;
  std::cout << "#adcs:         " << Signal.size()  << std::endl;
  std::cout << "Integral:      " << Integral      << std::endl; 
  std::cout << "PeakTime:      " << PeakTime      << std::endl;
  std::cout << "PeakAmplitude: " << PeakAmplitude << std::endl;
  std::cout << "Position:      " << "( " << Position.X() << ", " << Position.Y() << ", " << Position.Z() << ")" << std::endl;
  std::cout << "CNN:           " << "( " << CNN[0] << ", " << CNN[1] << ", " << CNN[2] << ")" << std::endl;
}

//********************************************************************
AnaTrajectoryPointPD::AnaTrajectoryPointPD(){
//********************************************************************

  Position_NoSCE  = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
  Direction_NoSCE = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
  Position        = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
  Direction       = TVector3(kFloatUnassigned,kFloatUnassigned,kFloatUnassigned);
}

//********************************************************************
AnaTrajectoryPointPD::AnaTrajectoryPointPD(const AnaTrajectoryPointPD& trp){
//********************************************************************

  Position_NoSCE  = trp.Position_NoSCE;
  Direction_NoSCE = trp.Direction_NoSCE;
  Position        = trp.Position;
  Direction       = trp.Direction;
}

//********************************************************************
void AnaTrajectoryPointPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrajectoryPointPD --------- " << std::endl;

  std::cout << "Position:      " << "( " << Position.X() << ", " << Position.Y() << ", " << Position.Z() << ")" << std::endl;
}

//********************************************************************
AnaParticlePD::AnaParticlePD():AnaParticle(){
//********************************************************************

  ParentID = kIntUnassigned;
  TrackID  = kIntUnassigned;
  ShowerID = kIntUnassigned;

  Type = kUnknown;
  IsCandidate = false;
  isBeamPart = false;
  isPandora  = false;
  BeamOrigin = false;

  FitPDG        = kFloatUnassigned;

  for(int i = 0; i < 4; i++){
    PositionStartSCE[i]=kFloatUnassigned;
    PositionEndSCE[i]=kFloatUnassigned;
  }
  for(int i = 0; i < 3; i++){
    DirectionStartSCE[i]=kFloatUnassigned;
    DirectionEndSCE[i]=kFloatUnassigned;
  }

  ThetaXZ = kFloatUnassigned;
  ThetaYZ = kFloatUnassigned;

  for (Int_t i=0;i<3;i++){
    NHitsPerPlane[i] = kIntUnassigned;
    truncLibo_dEdx=kFloatUnassigned;

    PIDA[i]=kFloatUnassigned;
    ReconPDG[i]=kIntUnassigned;
    for (int j=0; j<10; j++) {
      PID[i][j]=kFloatUnassigned;
      CALO[i][j]=kFloatUnassigned;
    }

    CNNscore[i]=kFloatUnassigned;
  }  
  vtx_CNN_michelscore=kFloatUnassigned;
  vtx_CNN_NHits=kIntUnassigned;

  Chi2Proton=kFloatUnassigned;
  Chi2Muon=kFloatUnassigned;
  Chi2ndf=kIntUnassigned;

  for (int i=0; i<2; i++){
    RangeMomentum[i] = kFloatUnassigned;
    RangeMomentum_alt[i] = kFloatUnassigned;    
  }
    
  Length_alt = kFloatUnassigned;
  Generation = kIntUnassigned;
  
  for (int i=0; i<3; i++){
    Hits[i].clear();
  }

  TrjPoints.clear();
}

//********************************************************************
AnaParticlePD::~AnaParticlePD(){
//********************************************************************

}

//********************************************************************
AnaParticlePD::AnaParticlePD(const AnaParticlePD& part):AnaParticle(part){
//********************************************************************

  ParentID       = part.ParentID;
  TrackID        = part.ParentID;
  ShowerID       = part.ParentID;

  Type           = part.Type;
  IsCandidate    = part.IsCandidate;
  isBeamPart     = part.isBeamPart;
  isPandora      = part.isPandora;
  BeamOrigin     = part.BeamOrigin;

  for(int i = 0; i < 4; i++)PositionStartSCE[i]=part.PositionStart[i];
  for(int i = 0; i < 3; i++)DirectionStartSCE[i]=part.DirectionStart[i];

  ThetaXZ = part.ThetaXZ;
  ThetaYZ = part.ThetaYZ;

  FitPDG         = part.FitPDG;

  for (Int_t i=0;i<3;i++){

    NHitsPerPlane[i] = part.NHitsPerPlane[i];

    PIDA[i]=part.PIDA[i];
    ReconPDG[i]=part.ReconPDG[i];
    
    for (int j=0; j<10; j++) {
      PID[i][j]=part.PID[i][j];
      CALO[i][j]=part.CALO[i][j];
    }

    CNNscore[i]=part.CNNscore[i];
  }
  vtx_CNN_michelscore=part.vtx_CNN_michelscore;
  vtx_CNN_NHits=part.vtx_CNN_NHits;

  Chi2Proton = part.Chi2Proton;
  Chi2Muon   = part.Chi2Muon;
  Chi2ndf    = part.Chi2ndf;

  truncLibo_dEdx = part.truncLibo_dEdx;
  
  for (int i=0; i<2; i++){
    RangeMomentum[i]     = part.RangeMomentum[i];  
    RangeMomentum_alt[i] = part.RangeMomentum_alt[i];  
  }


  Length_alt = part.Length_alt;
  Generation = part.Generation;

  for (int i=0; i<3; i++){
    Hits[i] = part.Hits[i];
  }

  TrjPoints = part.TrjPoints;
}

//********************************************************************
void AnaParticlePD::Print() const{
//********************************************************************

  AnaParticle::Print();
  
  std::cout << "-------- AnaParticlePD --------- " << std::endl;
  std::cout << "ParentID:                " << ParentID << std::endl;
  std::cout << "Type:                    " << Type << std::endl;
  std::cout << "Forced track ID:         " << TrackID << std::endl;
  std::cout << "Forced shower ID:        " << ShowerID << std::endl;
  std::cout << "isPandora:               " << isPandora << std::endl;
  std::cout << "BeamOrigin:              " << BeamOrigin << std::endl;
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

  std::cout << "CALO[0]:                 ";
  for (int i=0;i<10;i++) std::cout << CALO[0][i] << " ";
  std::cout << std::endl;

  std::cout << "CALO[1]:                 ";
  for (int i=0;i<10;i++) std::cout << CALO[1][i] << " ";
  std::cout << std::endl;


  std::cout << "CALO[2]:                 ";
  for (int i=0;i<10;i++) std::cout << CALO[2][i] << " ";
  std::cout << std::endl;

  std::cout << "Chi2Proton:              " << Chi2Proton << std::endl;
  std::cout << "Chi2Muon:                " << Chi2Muon << std::endl;
  std::cout << "Chi2ndf:                 " << Chi2ndf << std::endl;

  std::cout << "CNN score:               " << CNNscore[0] << " " << CNNscore[1] << " " << CNNscore[2] << std::endl;
  
  std::cout << "vtx CNN michelscore:     " << vtx_CNN_michelscore << std::endl;
  std::cout << "vtx CNN NHits:           " << vtx_CNN_NHits << std::endl;
  
  std::cout << "RangeMomentum:           " << RangeMomentum[0] << " " << RangeMomentum[1] << std::endl;

  std::cout << "Stored hits in plane 0,1,2:  " << Hits[0].size() << " " << Hits[1].size() << " " << Hits[2].size() << std::endl;    
}

//********************************************************************
AnaTrueParticlePD::AnaTrueParticlePD():AnaTrueParticle(){
//********************************************************************

  Pi0_decay_ID.clear();
  Origin  = kIntUnassigned;
  Generation  = kIntUnassigned;
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
  Generation = truePart.Generation;
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

//********************************************************************
AnaWireCNN::AnaWireCNN(){
//********************************************************************


  adcs.clear();
  wire=0;
  time=0;

}

//********************************************************************
AnaWireCNN::AnaWireCNN(const AnaWireCNN& CNNwire){
//********************************************************************

  adcs = CNNwire.adcs;
  wire = CNNwire.wire;
  time = CNNwire.time;
}

//********************************************************************
AnaBunchPD::AnaBunchPD():AnaBunch(){
//********************************************************************

  CNNwires.clear();
}

//********************************************************************
AnaBunchPD::~AnaBunchPD(){
//********************************************************************

}

//********************************************************************
AnaBunchPD::AnaBunchPD(const AnaBunchPD& bunch):AnaBunch(bunch){
//********************************************************************

  CNNwires = bunch.CNNwires;
}


//********************************************************************
void AnaBunchPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaBunchPD --------- " << std::endl;


  AnaBunch::Print();
  
  std::cout << "#wires:      " << CNNwires.size() << std::endl;  
}

//********************************************************************
AnaEventPD::AnaEventPD():AnaEvent(){
//********************************************************************

  CNNwires.clear();
}

//********************************************************************
AnaEventPD::~AnaEventPD(){
//********************************************************************

}

//********************************************************************
AnaEventPD::AnaEventPD(const AnaEventPD& event):AnaEvent(event){
//********************************************************************

  CNNwires = event.CNNwires;
}

//*****************************************************************************
AnaEventPD::AnaEventPD(const AnaSpillPD& spill, const AnaBunchPD& bunch):AnaEvent(spill,bunch) {
//*****************************************************************************

  //CNNwires = bunch.CNNwires;
}

//********************************************************************
void AnaEventPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaEventPD --------- " << std::endl;

  AnaEvent::Print();
  std::cout << "#wires:      " << CNNwires.size() << std::endl;
}

//********************************************************************
AnaEventInfoPD::AnaEventInfoPD():AnaEventInfo(){
//********************************************************************

  NominalBeamMom = kFloatUnassigned;
}

//********************************************************************
AnaEventInfoPD::~AnaEventInfoPD(){
//********************************************************************

}

//********************************************************************
AnaEventInfoPD::AnaEventInfoPD(const AnaEventInfoPD& eventinfo):AnaEventInfo(eventinfo){
//********************************************************************

  NominalBeamMom = eventinfo.NominalBeamMom;
}

//********************************************************************
void AnaEventInfoPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaEventInfoPD --------- " << std::endl;

  AnaEventInfo::Print();
  std::cout << "#NominalBeamMom:      " << NominalBeamMom << std::endl;
}
