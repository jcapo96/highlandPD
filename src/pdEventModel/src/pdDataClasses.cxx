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
  Position       = TVector3(0,0,0);
  Channel       = kIntUnassigned;
  View          = kIntUnassigned; 

  StartTick     = (UInt_t)kIntUnassigned;
  EndTick       = (UInt_t)kIntUnassigned; 
  
  dEdx          = kFloatUnassigned;
  dQdx          = kFloatUnassigned;
  dEdx_corr     = kFloatUnassigned;
  dQdx_corr     = kFloatUnassigned;
  ResidualRange = kFloatUnassigned;

  Signal.clear();
}

//********************************************************************
AnaHitPD::AnaHitPD(const AnaHitPD& hit){
//********************************************************************

  WireID.Plane  = hit.WireID.Plane;
  WireID.Wire   = hit.WireID.Wire; 
  Integral      = hit.Integral;     
  PeakTime      = hit.PeakTime;     
  PeakAmplitude = hit.PeakAmplitude;
  Position      = hit.Position;
  Channel       = hit.Channel;
  View          = hit.View;
  Signal        = hit.Signal;


  StartTick     = hit.StartTick;
  EndTick       = hit.EndTick  ; 
  
  dEdx           = hit.dEdx;
  dQdx           = hit.dQdx;
  dEdx_corr      = hit.dEdx_corr;
  dQdx_corr      = hit.dQdx_corr;
  ResidualRange  = hit.ResidualRange;
}

//********************************************************************
void AnaHitPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaHitPD --------- " << std::endl;

  std::cout << "WireID.Plane:  " << WireID.Plane  << std::endl;
  std::cout << "WireID.Wire:   " << WireID.Wire  << std::endl;
  std::cout << "View:          " << View << std::endl;
  std::cout << "#adcs:         " << Signal.size()  << std::endl;
  std::cout << "Integral:      " << Integral      << std::endl; 
  std::cout << "PeakTime:      " << PeakTime      << std::endl;
  std::cout << "PeakAmplitude: " << PeakAmplitude << std::endl;
  std::cout << "Position:      " << "( " << Position.X() << ", " << Position.Y() << ", " << Position.Z() << ")" << std::endl;
}

//********************************************************************
AnaParticlePD::AnaParticlePD():AnaParticle(){
//********************************************************************

  Type = kUnknown;
  isBeamPart    = false;
  isPandora = false;

  FitPDG        = kFloatUnassigned;

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

  Chi2Proton=kFloatUnassigned;
  Chi2ndf=kIntUnassigned;

  for (int i=0; i<2; i++){
    RangeMomentum[i] = kFloatUnassigned;
    RangeMomentum_alt[i] = kFloatUnassigned;    
  }
    
  Length_alt = kFloatUnassigned;

  for (int i=0; i<3; i++){
    Hits[i].clear();
  }
  
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
  isPandora      = part.isPandora;

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


  Chi2Proton = part.Chi2Proton;
  Chi2ndf    = part.Chi2ndf;

  truncLibo_dEdx = part.truncLibo_dEdx;
  
  for (int i=0; i<2; i++){
    RangeMomentum[i]     = part.RangeMomentum[i];  
    RangeMomentum_alt[i] = part.RangeMomentum_alt[i];  
  }


  Length_alt = part.Length_alt;
 
  for (int i=0; i<3; i++){
    Hits[i] = part.Hits[i];
  }

}

//********************************************************************
void AnaParticlePD::Print() const{
//********************************************************************

  AnaParticle::Print();
  
  std::cout << "-------- AnaParticlePD --------- " << std::endl;
  std::cout << "Type:                    " << Type << std::endl;
  std::cout << "isPandora:               " << isPandora << std::endl;
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
  std::cout << "Chi2ndf:                 " << Chi2ndf << std::endl;

  std::cout << "CNN score:               " << CNNscore[0] << " " << CNNscore[1] << " " << CNNscore[2] << std::endl;
  
  std::cout << "RangeMomentum:           " << RangeMomentum[0] << " " << RangeMomentum[1] << std::endl;

  std::cout << "Stored hits in plane 2:  " << Hits[2].size() << std::endl;    
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

//********************************************************************
AnaSpillPD::AnaSpillPD():AnaSpill(){
//********************************************************************

  ADC.clear();
  ADC.resize(15480);
  for (size_t w=0;w<ADC.size();w++)
    ADC[w].resize(6000,0);  



  reduced_adc_cnn_map.clear();
  reduced_adc_cnn_map_wires.clear();
  reduced_adc_cnn_map_times.clear();

}

//********************************************************************
AnaSpillPD::~AnaSpillPD(){
//********************************************************************

}

//********************************************************************
AnaSpillPD::AnaSpillPD(const AnaSpillPD& spill):AnaSpill(spill){
//********************************************************************

  ADC = spill.ADC;

  reduced_adc_cnn_map       = spill.reduced_adc_cnn_map      ;
  reduced_adc_cnn_map_wires = spill.reduced_adc_cnn_map_wires;
  reduced_adc_cnn_map_times = spill.reduced_adc_cnn_map_times;
}

//********************************************************************
void AnaSpillPD::RedoLinks(){
//********************************************************************

  for (size_t wi=0;wi<reduced_adc_cnn_map_wires.size();wi++){
    Int_t wire = reduced_adc_cnn_map_wires[wi];
    Int_t t0   = reduced_adc_cnn_map_times[wi];
    for (size_t s=0;s<reduced_adc_cnn_map[wi].size();s++){
      Int_t time = t0+s;
      ADC[wire][time]=reduced_adc_cnn_map[wi][s];
    }
  }

  // Redo the links in the bae class
  AnaSpill::RedoLinks();    
}

//********************************************************************
void AnaSpillPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaSpillPD --------- " << std::endl;

  AnaSpill::Print();

  std::cout << "ADC.size():            " << ADC.size() << std::endl;  
}

//********************************************************************
AnaEventPD::AnaEventPD():AnaEvent(){
//********************************************************************

  ADC.clear();
  ADC.resize(15480);
  for (size_t w=0;w<ADC.size();w++)
    ADC[w].resize(6000,0);
  
}

//********************************************************************
AnaEventPD::~AnaEventPD(){
//********************************************************************

}

//********************************************************************
AnaEventPD::AnaEventPD(const AnaEventPD& event):AnaEvent(event){
//********************************************************************

  
  ADC = event.ADC;
}

//*****************************************************************************
AnaEventPD::AnaEventPD(const AnaSpillPD& spill, const AnaBunch& bunch):AnaEvent(spill,bunch) {
//*****************************************************************************

  ADC.clear();
  ADC = spill.ADC;
}

//********************************************************************
void AnaEventPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaEventPD --------- " << std::endl;

  AnaEvent::Print();
  std::cout << "ADC.size():            " << ADC.size() << std::endl;  

}

