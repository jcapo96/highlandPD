#define pdDataClasses_C

#include "pdDataClasses.hxx"
#include "AnalysisUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdMomEstimation.hxx"
#include "pdMomReconstruction.hxx"
// #include "pdMomLikelihood.hxx" // Disabled for now, will implement later


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

  Distance_to_closest_particle = kFloatUnassigned;

  for (int i=0; i<3; i++){
    Hits[i].clear();
  }

  TrjPoints.clear();

  forced_daughter = false;
  forced_daughter_matched = false;
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

  Distance_to_closest_particle = part.Distance_to_closest_particle;

  for (int i=0; i<3; i++){
    Hits[i] = part.Hits[i];
  }

  TrjPoints = part.TrjPoints;

  forced_daughter = part.forced_daughter;
  forced_daughter_matched = part.forced_daughter_matched;
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
AnaVertexFittedParticlePD::AnaVertexFittedParticlePD():AnaParticlePD(){
//********************************************************************

  for(int i = 0; i < 3; i++){
    FitPosition[i] = kFloatUnassigned;
    FitDirection[i] = kFloatUnassigned;
  }
}

//********************************************************************
AnaVertexFittedParticlePD::~AnaVertexFittedParticlePD(){
//********************************************************************

}

//********************************************************************
AnaVertexFittedParticlePD::AnaVertexFittedParticlePD(const AnaVertexFittedParticlePD& fittedPart):AnaParticlePD(fittedPart){
//********************************************************************

  for(int i = 0; i < 3; i++){
    FitPosition[i] = fittedPart.FitPosition[i];
    FitDirection[i] = fittedPart.FitDirection[i];
  }
}

//********************************************************************
void AnaVertexFittedParticlePD::Print() const{
//********************************************************************

  AnaParticlePD::Print();

  std::cout << "-------- AnaVertexFittedParticlePD --------- " << std::endl;
  std::cout << "FitPosition: (" << FitPosition[0] << ", " << FitPosition[1] << ", " << FitPosition[2] << ")" << std::endl;
  std::cout << "FitDirection: (" << FitDirection[0] << ", " << FitDirection[1] << ", " << FitDirection[2] << ")" << std::endl;
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
  EmptyEvent = false;
  HasPandora = false;
}

//********************************************************************
AnaEventInfoPD::~AnaEventInfoPD(){
//********************************************************************

}

//********************************************************************
AnaEventInfoPD::AnaEventInfoPD(const AnaEventInfoPD& eventinfo):AnaEventInfo(eventinfo){
//********************************************************************

  NominalBeamMom = eventinfo.NominalBeamMom;
  EmptyEvent = eventinfo.EmptyEvent;
  HasPandora = eventinfo.HasPandora;
}

//********************************************************************
void AnaEventInfoPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaEventInfoPD --------- " << std::endl;

  AnaEventInfo::Print();
  std::cout << "#NominalBeamMom:      " << NominalBeamMom << std::endl;
}

//********************************************************************
AnaTrueEquivalentVertexPD::AnaTrueEquivalentVertexPD(){
//********************************************************************

  TrueParticles.clear();
  OriginalDistance = kFloatUnassigned;
  MinimumDistance = kFloatUnassigned;
  OpeningAngle = kFloatUnassigned;
  Position[0] = kFloatUnassigned;
  Position[1] = kFloatUnassigned;
  Position[2] = kFloatUnassigned;
  Direction[0] = kFloatUnassigned;
  Direction[1] = kFloatUnassigned;
  Direction[2] = kFloatUnassigned;
  FitParticles.clear();
  FitDirection[0] = kFloatUnassigned;
  FitDirection[1] = kFloatUnassigned;
  FitDirection[2] = kFloatUnassigned;
  PositionPandora[0] = kFloatUnassigned;
  PositionPandora[1] = kFloatUnassigned;
  PositionPandora[2] = kFloatUnassigned;
  DegeneracyBeforeScoring = 0;
  DegeneracyAfterScoring = 0;
  NRecoParticles = 0;
}

//********************************************************************
AnaTrueEquivalentVertexPD::~AnaTrueEquivalentVertexPD(){
//********************************************************************

}

//********************************************************************
AnaTrueEquivalentVertexPD::AnaTrueEquivalentVertexPD(const AnaTrueEquivalentVertexPD& vertex){
//********************************************************************

  TrueParticles = vertex.TrueParticles;
  OriginalDistance = vertex.OriginalDistance;
  MinimumDistance = vertex.MinimumDistance;
  OpeningAngle = vertex.OpeningAngle;
  Position[0] = vertex.Position[0];
  Position[1] = vertex.Position[1];
  Position[2] = vertex.Position[2];
  Direction[0] = vertex.Direction[0];
  Direction[1] = vertex.Direction[1];
  Direction[2] = vertex.Direction[2];
  FitParticles = vertex.FitParticles;
  FitDirection[0] = vertex.FitDirection[0];
  FitDirection[1] = vertex.FitDirection[1];
  FitDirection[2] = vertex.FitDirection[2];
  PositionPandora[0] = vertex.PositionPandora[0];
  PositionPandora[1] = vertex.PositionPandora[1];
  PositionPandora[2] = vertex.PositionPandora[2];
  DegeneracyBeforeScoring = vertex.DegeneracyBeforeScoring;
  DegeneracyAfterScoring = vertex.DegeneracyAfterScoring;
  NRecoParticles = vertex.NRecoParticles;
}

//********************************************************************
void AnaTrueEquivalentVertexPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrueEquivalentVertexPD --------- " << std::endl;
  std::cout << "TrueParticles size:    " << TrueParticles.size() << std::endl;
  std::cout << "OriginalDistance:      " << OriginalDistance << " cm" << std::endl;
  std::cout << "MinimumDistance:       " << MinimumDistance << " cm" << std::endl;
  std::cout << "OpeningAngle:          " << OpeningAngle << " degrees" << std::endl;
  std::cout << "Position:              " << Position[0] << " " << Position[1] << " " << Position[2] << std::endl;
  std::cout << "Direction:             " << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
}

//********************************************************************
AnaVertexPD::AnaVertexPD():AnaVertexB(){
//********************************************************************

  UniqueID = kIntUnassigned;
  NParticles = 0;
  Particles.clear();
  OriginalDistance = kFloatUnassigned;
  Position[0] = kFloatUnassigned;
  Position[1] = kFloatUnassigned;
  Position[2] = kFloatUnassigned;
  Momentum[0] = kFloatUnassigned;
  Momentum[1] = kFloatUnassigned;
  Momentum[2] = kFloatUnassigned;
  Direction[0] = kFloatUnassigned;
  Direction[1] = kFloatUnassigned;
  Direction[2] = kFloatUnassigned;
  Energy = kFloatUnassigned;
  OpeningAngle = kFloatUnassigned;
  AngleWithBeam = kFloatUnassigned;
  NPotentialParents = kIntUnassigned;
  Generation = kIntUnassigned;
  Process = kIntUnassigned;
  MinimumDistance = kFloatUnassigned;
  Score = kFloatUnassigned;
  ParentID = kIntUnassigned;
  FitParticles.clear();
  FitDirection[0] = kFloatUnassigned;
  FitDirection[1] = kFloatUnassigned;
  FitDirection[2] = kFloatUnassigned;
  PositionPandora[0] = kFloatUnassigned;
  PositionPandora[1] = kFloatUnassigned;
  PositionPandora[2] = kFloatUnassigned;
  DegeneracyBeforeScoring = 0;
  DegeneracyAfterScoring = 0;
  NRecoParticles = 0;
}

//********************************************************************
AnaVertexPD::~AnaVertexPD(){
//********************************************************************

}

//********************************************************************
AnaVertexPD::AnaVertexPD(const AnaVertexPD& vertex):AnaVertexB(vertex){
//********************************************************************

  UniqueID = vertex.UniqueID;
  NParticles = vertex.NParticles;
  Particles = vertex.Particles;
  OriginalDistance = vertex.OriginalDistance;
  Position[0] = vertex.Position[0];
  Position[1] = vertex.Position[1];
  Position[2] = vertex.Position[2];
  Momentum[0] = vertex.Momentum[0];
  Momentum[1] = vertex.Momentum[1];
  Momentum[2] = vertex.Momentum[2];
  Direction[0] = vertex.Direction[0];
  Direction[1] = vertex.Direction[1];
  Direction[2] = vertex.Direction[2];
  Energy = vertex.Energy;
  OpeningAngle = vertex.OpeningAngle;
  AngleWithBeam = vertex.AngleWithBeam;
  NPotentialParents = vertex.NPotentialParents;
  Generation = vertex.Generation;
  Process = vertex.Process;
  MinimumDistance = vertex.MinimumDistance;
  Score = vertex.Score;
  ParentID = vertex.ParentID;
  FitParticles = vertex.FitParticles;
  FitDirection[0] = vertex.FitDirection[0];
  FitDirection[1] = vertex.FitDirection[1];
  FitDirection[2] = vertex.FitDirection[2];
  PositionPandora[0] = vertex.PositionPandora[0];
  PositionPandora[1] = vertex.PositionPandora[1];
  PositionPandora[2] = vertex.PositionPandora[2];
  DegeneracyBeforeScoring = vertex.DegeneracyBeforeScoring;
  DegeneracyAfterScoring = vertex.DegeneracyAfterScoring;
  NRecoParticles = vertex.NRecoParticles;
}

//********************************************************************
void AnaVertexPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaVertexPD --------- " << std::endl;

  AnaVertexB::Print();
  std::cout << "UniqueID:              " << UniqueID << std::endl;
  std::cout << "NParticles:            " << NParticles << std::endl;
  std::cout << "Particles size:        " << Particles.size() << std::endl;
  std::cout << "OriginalDistance:      " << OriginalDistance << " cm" << std::endl;
  std::cout << "Position:              " << Position[0] << " " << Position[1] << " " << Position[2] << std::endl;
  std::cout << "Momentum:              " << Momentum[0] << " " << Momentum[1] << " " << Momentum[2] << std::endl;
  std::cout << "Direction:             " << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
  std::cout << "Energy:                " << Energy << " GeV" << std::endl;
  std::cout << "OpeningAngle:          " << OpeningAngle << " rad" << std::endl;
  std::cout << "AngleWithBeam:         " << AngleWithBeam << " rad" << std::endl;
  std::cout << "NPotentialParents:      " << NPotentialParents << std::endl;
  std::cout << "Generation:            " << Generation << std::endl;
  std::cout << "Process:               " << Process << std::endl;
  std::cout << "MinimumDistance:       " << MinimumDistance << " cm" << std::endl;
  std::cout << "Score:                 " << Score << std::endl;
  std::cout << "ParentID:              " << ParentID << std::endl;
}

//********************************************************************
void AnaVertexPD::EnsureParticleMomentum(){

  // Check each particle in the vertex
  for (size_t i = 0; i < Particles.size(); i++) {
    AnaParticlePD* particle = Particles[i];
    if (!particle) continue;

    // Check if particle already has valid momentum
    if (particle->Momentum > 0 && particle->Momentum != -999) {
      // Particle already has valid momentum, skip
      continue;
    }

    Float_t calculatedMomentum = -999;

    // Priority 1: Calorimetric method (best for interacting pions from K0 decay)
    // Sums all deposited energy from particle + daughters using pitch-corrected path lengths
    if (!particle->Hits[2].empty()) {
      calculatedMomentum = pdMomReconstruction::EstimateMomentumCalorimetric(particle, 211);
    }

    // Priority 2: Track-length extension method (fallback for through-going pions)
    // This method fits the dE/dx shape to estimate how much further the particle would travel
    if (calculatedMomentum <= 0 || calculatedMomentum == -999) {
      if (!particle->Hits[2].empty() && particle->Length > 0) {
        calculatedMomentum = pdMomEstimation::EstimateMomentumWithExtension(particle, 211);
      }
    }

    // Fallback methods if calorimetric and extension methods both fail
    if (calculatedMomentum <= 0 || calculatedMomentum == -999) {

      // Calculate RangeMomentum if not available from input tree
      if ((particle->RangeMomentum[0] == -999 || particle->RangeMomentum[0] <= 0) &&
          particle->Length > 0 && particle->Length != -999) {
        particle->RangeMomentum[0] = pdAnaUtils::ComputeRangeMomentum(particle->Length, 2212);
      }
      if ((particle->RangeMomentum[1] == -999 || particle->RangeMomentum[1] <= 0) &&
          particle->Length > 0 && particle->Length != -999) {
        particle->RangeMomentum[1] = pdAnaUtils::ComputeRangeMomentum(particle->Length, 13);
      }

      // Priority 2: Use average of proton and muon range momenta as approximation for pion
      if (particle->RangeMomentum[0] > 0 && particle->RangeMomentum[0] != -999 &&
          particle->RangeMomentum[1] > 0 && particle->RangeMomentum[1] != -999) {
        calculatedMomentum = (particle->RangeMomentum[0] + particle->RangeMomentum[1]) / 2.0;
      }
      // Priority 3: Use muon range momentum if available (closer to pion mass than proton)
      else if (particle->RangeMomentum[1] > 0 && particle->RangeMomentum[1] != -999) {
        calculatedMomentum = particle->RangeMomentum[1];
      }
      // Priority 4: Use proton range momentum if available
      else if (particle->RangeMomentum[0] > 0 && particle->RangeMomentum[0] != -999) {
        calculatedMomentum = particle->RangeMomentum[0];
      }
    }

    // If we successfully calculated momentum, assign it
    if (calculatedMomentum > 0 && calculatedMomentum != -999) {
      particle->Momentum = calculatedMomentum;
    }

  }
}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(): AnaParticleB(){
//********************************************************************

  UniqueID = kIntUnassigned;
  Vertex = NULL;
  Parent = NULL;
  TrueEquivalentNeutralParticle = NULL;
  ImpactParameter = kFloatUnassigned;
  Mass = kFloatUnassigned;
  Momentum = kFloatUnassigned;
  PDG = kIntUnassigned;
  Lifetime = kFloatUnassigned;
  DecayLength = kFloatUnassigned;
  NRecoHitsInVertex = kIntUnassigned;
  NeutralScore = kFloatUnassigned;
  HitsAlignment = kFloatUnassigned;
  NHitsInCylinder = kIntUnassigned;
  HitsAvgDistance = kFloatUnassigned;
  HitsRMSDistance = kFloatUnassigned;
  HitsLongitudinalSpan = kFloatUnassigned;
  FitParent = NULL;
}

//********************************************************************
AnaNeutralParticlePD::~AnaNeutralParticlePD(){
//********************************************************************

}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(const AnaNeutralParticlePD& neutralParticle): AnaParticleB(neutralParticle){
//********************************************************************

  UniqueID = neutralParticle.UniqueID;
  Vertex = neutralParticle.Vertex;
  Parent = neutralParticle.Parent;
  TrueEquivalentNeutralParticle = neutralParticle.TrueEquivalentNeutralParticle;
  RecoParticle = neutralParticle.RecoParticle;
  ImpactParameter = neutralParticle.ImpactParameter;
  Mass = neutralParticle.Mass;
  Momentum = neutralParticle.Momentum;
  PDG = neutralParticle.PDG;
  Lifetime = neutralParticle.Lifetime;
  DecayLength = neutralParticle.DecayLength;
  NRecoHitsInVertex = neutralParticle.NRecoHitsInVertex;
  NeutralScore = neutralParticle.NeutralScore;
  HitsAlignment = neutralParticle.HitsAlignment;
  NHitsInCylinder = neutralParticle.NHitsInCylinder;
  HitsAvgDistance = neutralParticle.HitsAvgDistance;
  HitsRMSDistance = neutralParticle.HitsRMSDistance;
  HitsLongitudinalSpan = neutralParticle.HitsLongitudinalSpan;
  FitParent = neutralParticle.FitParent;
}

//********************************************************************
void AnaNeutralParticlePD::Print() const{
//********************************************************************

  std::cout << "-------- AnaNeutralParticlePD --------- " << std::endl;

  AnaParticleB::Print();

  std::cout << "UniqueID:              " << UniqueID << std::endl;
  std::cout << "Vertex:                " << (Vertex ? "Yes" : "No") << std::endl;
  std::cout << "Parent:                " << (Parent ? "Yes" : "No") << std::endl;
  std::cout << "TrueEquivalentNeutralParticle:   " << (TrueEquivalentNeutralParticle ? "Yes" : "No") << std::endl;
  std::cout << "RecoParticle:          " << (RecoParticle ? "Yes" : "No") << std::endl;
  std::cout << "ImpactParameter:       " << ImpactParameter << " cm" << std::endl;
  std::cout << "Mass:                  " << Mass << " GeV/c²" << std::endl;
  std::cout << "Momentum:              " << Momentum << " GeV/c" << std::endl;
  std::cout << "PDG:                   " << PDG << std::endl;
  std::cout << "Lifetime:              " << Lifetime << " ns" << std::endl;
  std::cout << "DecayLength:           " << DecayLength << " cm" << std::endl;
  std::cout << "NRecoHitsInVertex:     " << NRecoHitsInVertex << std::endl;
  std::cout << "RecoParticle:          " << (RecoParticle ? "Yes" : "No") << std::endl;
}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::AnaTrueEquivalentNeutralParticlePD(){
//********************************************************************

  TrueEquivalentVertex = NULL;
  TrueParent = NULL;
  Position[0] = kFloatUnassigned;
  Position[1] = kFloatUnassigned;
  Position[2] = kFloatUnassigned;
  Direction[0] = kFloatUnassigned;
  Direction[1] = kFloatUnassigned;
  Direction[2] = kFloatUnassigned;
  PositionEnd[0] = kFloatUnassigned;
  PositionEnd[1] = kFloatUnassigned;
  PositionEnd[2] = kFloatUnassigned;
  DirectionEnd[0] = kFloatUnassigned;
  DirectionEnd[1] = kFloatUnassigned;
  DirectionEnd[2] = kFloatUnassigned;
  Length = kFloatUnassigned;
  Momentum = kFloatUnassigned;
  PDG = kIntUnassigned;
  Generation = kIntUnassigned;
  Process = kIntUnassigned;
  Mass = kFloatUnassigned;
}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::~AnaTrueEquivalentNeutralParticlePD(){
//********************************************************************

}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::AnaTrueEquivalentNeutralParticlePD(const AnaTrueEquivalentNeutralParticlePD& trueEquivalentNeutralPart){
//********************************************************************

  TrueEquivalentVertex = trueEquivalentNeutralPart.TrueEquivalentVertex;
  TrueParent = trueEquivalentNeutralPart.TrueParent;
  Position[0] = trueEquivalentNeutralPart.Position[0];
  Position[1] = trueEquivalentNeutralPart.Position[1];
  Position[2] = trueEquivalentNeutralPart.Position[2];
  Direction[0] = trueEquivalentNeutralPart.Direction[0];
  Direction[1] = trueEquivalentNeutralPart.Direction[1];
  Direction[2] = trueEquivalentNeutralPart.Direction[2];
  PositionEnd[0] = trueEquivalentNeutralPart.PositionEnd[0];
  PositionEnd[1] = trueEquivalentNeutralPart.PositionEnd[1];
  PositionEnd[2] = trueEquivalentNeutralPart.PositionEnd[2];
  DirectionEnd[0] = trueEquivalentNeutralPart.DirectionEnd[0];
  DirectionEnd[1] = trueEquivalentNeutralPart.DirectionEnd[1];
  DirectionEnd[2] = trueEquivalentNeutralPart.DirectionEnd[2];
  Length = trueEquivalentNeutralPart.Length;
  Momentum = trueEquivalentNeutralPart.Momentum;
  PDG = trueEquivalentNeutralPart.PDG;
  Generation = trueEquivalentNeutralPart.Generation;
  Process = trueEquivalentNeutralPart.Process;
  Mass = trueEquivalentNeutralPart.Mass;
}

//********************************************************************
void AnaTrueEquivalentNeutralParticlePD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrueNeutralParticlePD --------- " << std::endl;
  std::cout << "TrueEquivalentVertex:   " << (TrueEquivalentVertex ? "Set" : "NULL") << std::endl;
  std::cout << "TrueParent:             " << (TrueParent ? "Set" : "NULL") << std::endl;
  std::cout << "Position:               " << Position[0] << " " << Position[1] << " " << Position[2] << std::endl;
  std::cout << "Direction:              " << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
  std::cout << "PositionEnd:            " << PositionEnd[0] << " " << PositionEnd[1] << " " << PositionEnd[2] << std::endl;
  std::cout << "DirectionEnd:           " << DirectionEnd[0] << " " << DirectionEnd[1] << " " << DirectionEnd[2] << std::endl;
  std::cout << "Length:                 " << Length << " cm" << std::endl;
  std::cout << "Momentum:               " << Momentum << " GeV/c" << std::endl;
  std::cout << "PDG:                    " << PDG << std::endl;
  std::cout << "Generation:             " << Generation << std::endl;
  std::cout << "Process:                " << Process << std::endl;
  std::cout << "Mass:                   " << Mass << " GeV/c²" << std::endl;
}

