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
  TrajectoryDirection = TVector3(kFloatUnassigned, kFloatUnassigned, kFloatUnassigned);
  TrajectoryDirectionNPoints = 0;

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
  TrajectoryDirection = part.TrajectoryDirection;
  TrajectoryDirectionNPoints = part.TrajectoryDirectionNPoints;

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
  MinimumDistancePandora = kFloatUnassigned;
  MinimumDistanceFit = kFloatUnassigned;
  ClosestPoint1[0] = ClosestPoint1[1] = ClosestPoint1[2] = kFloatUnassigned;
  ClosestPoint2[0] = ClosestPoint2[1] = ClosestPoint2[2] = kFloatUnassigned;
  ClosestPoint1Pandora[0] = ClosestPoint1Pandora[1] = ClosestPoint1Pandora[2] = kFloatUnassigned;
  ClosestPoint2Pandora[0] = ClosestPoint2Pandora[1] = ClosestPoint2Pandora[2] = kFloatUnassigned;
  ClosestPoint1Fit[0] = ClosestPoint1Fit[1] = ClosestPoint1Fit[2] = kFloatUnassigned;
  ClosestPoint2Fit[0] = ClosestPoint2Fit[1] = ClosestPoint2Fit[2] = kFloatUnassigned;
  Score = kFloatUnassigned;
  ScorePandora = kFloatUnassigned;
  ScoreFit = kFloatUnassigned;
  ParentID = kIntUnassigned;
  DirectionFit[0] = kFloatUnassigned;
  DirectionFit[1] = kFloatUnassigned;
  DirectionFit[2] = kFloatUnassigned;
  PositionPandora[0] = kFloatUnassigned;
  PositionPandora[1] = kFloatUnassigned;
  PositionPandora[2] = kFloatUnassigned;
  PositionFit[0] = kFloatUnassigned;
  PositionFit[1] = kFloatUnassigned;
  PositionFit[2] = kFloatUnassigned;
  DirectionFit[0] = kFloatUnassigned;
  DirectionFit[1] = kFloatUnassigned;
  DirectionFit[2] = kFloatUnassigned;
  IsJustAverage = 0;
  Degeneracy = 0;
  for (Int_t i = 0; i < 5; ++i) {
    DegDist[i] = kFloatUnassigned;
  }
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
  MinimumDistancePandora = vertex.MinimumDistancePandora;
  MinimumDistanceFit = vertex.MinimumDistanceFit;
  ClosestPoint1[0] = vertex.ClosestPoint1[0];
  ClosestPoint1[1] = vertex.ClosestPoint1[1];
  ClosestPoint1[2] = vertex.ClosestPoint1[2];
  ClosestPoint2[0] = vertex.ClosestPoint2[0];
  ClosestPoint2[1] = vertex.ClosestPoint2[1];
  ClosestPoint2[2] = vertex.ClosestPoint2[2];
  ClosestPoint1Pandora[0] = vertex.ClosestPoint1Pandora[0];
  ClosestPoint1Pandora[1] = vertex.ClosestPoint1Pandora[1];
  ClosestPoint1Pandora[2] = vertex.ClosestPoint1Pandora[2];
  ClosestPoint2Pandora[0] = vertex.ClosestPoint2Pandora[0];
  ClosestPoint2Pandora[1] = vertex.ClosestPoint2Pandora[1];
  ClosestPoint2Pandora[2] = vertex.ClosestPoint2Pandora[2];
  ClosestPoint1Fit[0] = vertex.ClosestPoint1Fit[0];
  ClosestPoint1Fit[1] = vertex.ClosestPoint1Fit[1];
  ClosestPoint1Fit[2] = vertex.ClosestPoint1Fit[2];
  ClosestPoint2Fit[0] = vertex.ClosestPoint2Fit[0];
  ClosestPoint2Fit[1] = vertex.ClosestPoint2Fit[1];
  ClosestPoint2Fit[2] = vertex.ClosestPoint2Fit[2];
  Score = vertex.Score;
  ScorePandora = vertex.ScorePandora;
  ScoreFit = vertex.ScoreFit;
  ParentID = vertex.ParentID;
  DirectionFit[0] = vertex.DirectionFit[0];
  DirectionFit[1] = vertex.DirectionFit[1];
  DirectionFit[2] = vertex.DirectionFit[2];
  PositionPandora[0] = vertex.PositionPandora[0];
  PositionPandora[1] = vertex.PositionPandora[1];
  PositionPandora[2] = vertex.PositionPandora[2];
  PositionFit[0] = vertex.PositionFit[0];
  PositionFit[1] = vertex.PositionFit[1];
  PositionFit[2] = vertex.PositionFit[2];
  DirectionFit[0] = vertex.DirectionFit[0];
  DirectionFit[1] = vertex.DirectionFit[1];
  DirectionFit[2] = vertex.DirectionFit[2];
  IsJustAverage = vertex.IsJustAverage;
  Degeneracy = vertex.Degeneracy;
  for (Int_t i = 0; i < 5; ++i) {
    DegDist[i] = vertex.DegDist[i];
  }
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
AnaCreationVertexPD::AnaCreationVertexPD() : AnaVertexPD(){
//********************************************************************
  BeamParticle = NULL;
  SecondParticle = NULL;
  ProtonScore = -999.0;
  DistanceScore = -999.0;
  MinDistanceScore = -999.0;
  Position[0] = -999.0;
  Position[1] = -999.0;
  Position[2] = -999.0;
  ClosestPointBeam[0] = kFloatUnassigned;
  ClosestPointBeam[1] = kFloatUnassigned;
  ClosestPointBeam[2] = kFloatUnassigned;
  ClosestPointSecond[0] = kFloatUnassigned;
  ClosestPointSecond[1] = kFloatUnassigned;
  ClosestPointSecond[2] = kFloatUnassigned;
}

//********************************************************************
AnaCreationVertexPD::~AnaCreationVertexPD(){
//********************************************************************
}

//********************************************************************
AnaAnnihilationVertexPD::AnaAnnihilationVertexPD(){
//********************************************************************
  UniqueID = kIntUnassigned;
  NParticles = 0;
  Particles.clear();
  PositionPandora[0] = kFloatUnassigned;
  PositionPandora[1] = kFloatUnassigned;
  PositionPandora[2] = kFloatUnassigned;
  PositionFit[0] = kFloatUnassigned;
  PositionFit[1] = kFloatUnassigned;
  PositionFit[2] = kFloatUnassigned;
  ClosestPointPandora1[0] = kFloatUnassigned;
  ClosestPointPandora1[1] = kFloatUnassigned;
  ClosestPointPandora1[2] = kFloatUnassigned;
  ClosestPointPandora2[0] = kFloatUnassigned;
  ClosestPointPandora2[1] = kFloatUnassigned;
  ClosestPointPandora2[2] = kFloatUnassigned;
  ClosestPointFit1[0] = kFloatUnassigned;
  ClosestPointFit1[1] = kFloatUnassigned;
  ClosestPointFit1[2] = kFloatUnassigned;
  ClosestPointFit2[0] = kFloatUnassigned;
  ClosestPointFit2[1] = kFloatUnassigned;
  ClosestPointFit2[2] = kFloatUnassigned;
  MinimumDistancePandora = kFloatUnassigned;
  MinimumDistanceFit = kFloatUnassigned;
  OriginalDistance = kFloatUnassigned;
  Degeneracy = kIntUnassigned;
  Momentum = kFloatUnassigned;
  InvariantMass = kFloatUnassigned;
  MomentumPandora = kFloatUnassigned;
  InvariantMassPandora = kFloatUnassigned;
  MomentumFit = kFloatUnassigned;
  InvariantMassFit = kFloatUnassigned;
}

//********************************************************************
AnaAnnihilationVertexPD::~AnaAnnihilationVertexPD(){
//********************************************************************
}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(): AnaParticleB(){
//********************************************************************

  UniqueID = kIntUnassigned;
  AnnihilationVertex = NULL;
  CreationVertex = NULL;
  Parent = NULL;
  LengthPandora = kFloatUnassigned;
  LengthFit = kFloatUnassigned;
  AlignmentPandora = kFloatUnassigned;
  AlignmentFit = kFloatUnassigned;
}

//********************************************************************
AnaNeutralParticlePD::~AnaNeutralParticlePD(){
//********************************************************************

}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(const AnaNeutralParticlePD& neutralParticle): AnaParticleB(neutralParticle){
//********************************************************************

  UniqueID = neutralParticle.UniqueID;
  AnnihilationVertex = neutralParticle.AnnihilationVertex;
  CreationVertex = neutralParticle.CreationVertex;
  Parent = neutralParticle.Parent;
  LengthPandora = neutralParticle.LengthPandora;
  LengthFit = neutralParticle.LengthFit;
  AlignmentPandora = neutralParticle.AlignmentPandora;
  AlignmentFit = neutralParticle.AlignmentFit;
}

//********************************************************************
void AnaNeutralParticlePD::Print() const{
//********************************************************************

  std::cout << "-------- AnaNeutralParticlePD --------- " << std::endl;

  AnaParticleB::Print();

  std::cout << "UniqueID:              " << UniqueID << std::endl;
  std::cout << "AnnihilationVertex:    " << (AnnihilationVertex ? "Yes" : "No") << std::endl;
  std::cout << "CreationVertex:        " << (CreationVertex ? "Yes" : "No") << std::endl;
  std::cout << "Parent:                " << (Parent ? "Yes" : "No") << std::endl;
  std::cout << "LengthPandora:         " << LengthPandora << " cm" << std::endl;
  std::cout << "LengthFit:             " << LengthFit << " cm" << std::endl;
  std::cout << "AlignmentPandora:      " << AlignmentPandora << " rad" << std::endl;
  std::cout << "AlignmentFit:          " << AlignmentFit << " rad" << std::endl;
}

