#include "CNNUtils.hxx"
#include "timeUtils.hxx"
#include "Parameters.hxx"
#include <chrono>
#include <thread>


#ifdef  CompileTF
#include "tensorflow/core/public/session.h"  // anselmo
#endif


std::string wspacesA(size_t n){
  std::string ws="";
  for (size_t i=0;i<n;i++)
    ws +=" ";
  return ws;
}


/*
  produce
    setWireDriftDataFromHit;    
      resizeView
      setWireData
    predictIdVectors;
      bufferPatch
        patchFromOriginalView
          downscale
            downscaleMean
              scaleAdcSamples
                scaleAdcSample
      fNNet->Run(inps);

*/



std::pair<Int_t, Int_t> wire_mapping[12][3]=    {{std::pair<Int_t,Int_t>(0,799),   
                                                  std::pair<Int_t,Int_t>(800  , 1599),   
                                                  std::pair<Int_t,Int_t>(1600 , 2079)},   
                                                 {std::pair<Int_t,Int_t>(0    , 799),    
                                                  std::pair<Int_t,Int_t>(800  , 1599),   
                                                  std::pair<Int_t,Int_t>(2080 , 2559)},   
                                                 {std::pair<Int_t,Int_t>(2560 , 3359),   
                                                  std::pair<Int_t,Int_t>(3360 , 4159),   
                                                  std::pair<Int_t,Int_t>(4160 , 4639)},   
                                                 {std::pair<Int_t,Int_t>(2560 , 3359),   
                                                  std::pair<Int_t,Int_t>(3360 , 4159),   
                                                  std::pair<Int_t,Int_t>(4640 , 5119)},   
                                                 {std::pair<Int_t,Int_t>(5120 , 5919),   
                                                  std::pair<Int_t,Int_t>(5920 , 6719),   
                                                  std::pair<Int_t,Int_t>(6720 , 7199)},   
                                                 {std::pair<Int_t,Int_t>(5120 , 5919),   
                                                  std::pair<Int_t,Int_t>(5920 , 6719),   
                                                  std::pair<Int_t,Int_t>(7200 , 7679)},   
                                                 {std::pair<Int_t,Int_t>(7680 , 8479),   
                                                  std::pair<Int_t,Int_t>(8480 , 9279),   
                                                  std::pair<Int_t,Int_t>(9280 , 9759)},   
                                                 {std::pair<Int_t,Int_t>(7680 , 8479),   
                                                  std::pair<Int_t,Int_t>(8480 , 9279),   
                                                  std::pair<Int_t,Int_t>(9760 , 10239)}, 
                                                 {std::pair<Int_t,Int_t>(10240 , 11039),
                                                  std::pair<Int_t,Int_t>(11040 , 11839),
                                                  std::pair<Int_t,Int_t>(11840 , 12319)},
                                                 {std::pair<Int_t,Int_t>(10240 , 11039),
                                                  std::pair<Int_t,Int_t>(11040 , 11839),
                                                  std::pair<Int_t,Int_t>(12320 , 12799)},
                                                 {std::pair<Int_t,Int_t>(12800 , 13599),
                                                  std::pair<Int_t,Int_t>(13600 , 14399),
                                                  std::pair<Int_t,Int_t>(14400 , 14879)},
                                                 {std::pair<Int_t,Int_t>(12800 , 13599),
                                                  std::pair<Int_t,Int_t>(13600 , 14399),
                                                  std::pair<Int_t,Int_t>(14880 , 15359)}};


//*******************************************************
CNNUtils::CNNUtils():
  fDownscaleMode(kMax),
  fDriftWindow(6), 
  fDownscaleFullView(false),
  // set fixed threshold of 10 ADC counts for counting the sum

    //  , fNoiseSigma(0)
    //  , fCoherentSigma(0)
  fCalibrateAmpl(true),
  fCalibrateLifetime(false),
  fAdcSumOverThr(0),
  fAdcSumThr(10),
  fAdcAreaOverThr(0),
  fPatchSizeW(48), 
  fPatchSizeD(48),
  fCurrentWireIdx(99999),
  fCurrentScaledDrift(99999),
  fNNet(0)
{
  /*

  size_t fDriftWindow;
  bool fDownscaleFullView;
  std::vector<float> fAmplCalibConst;
  bool fCalibrateAmpl;
  bool fCalibrateLifetime;
  float fAdcMax, fAdcMin, fAdcScale, fAdcOffset, fAdcZero;
  double fAdcSumOverThr;
  double fAdcSumThr;
  size_t fAdcAreaOverThr;
  size_t fPatchSizeW, fPatchSizeD;
  mutable size_t fCurrentWireIdx, fCurrentScaledDrift;
  */
  
//*******************************************************


  fNNetModelFilePath =   (std::string)getenv("PDUTILSROOT")+"/data/"+ND::params().GetParameterS("pdUtils.CNN.PBfile");
  //  fNNetOutputs = "";
  
  if (fNNet) delete fNNet;
  fNNet = 0;


  if ((fNNetModelFilePath.length() > 3) &&
      (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 3, 3, ".pb") == 0)) {    
    fNNet = new nnet2::TfModelInterface(fNNetModelFilePath.c_str());
  }
  else {
    std::cout  << "File name extension not supported. Only .pb files are supported" << std::endl;
  }

  if (!fNNet) { std:: cout << "Loading model from file failed." << std::endl;exit(1);}

    
  //    fAmplCalibConst.resize(fGeometry->MaxPlanes());
  fAmplCalibConst.resize(3);
  if (fCalibrateAmpl) {
    std::cout << "Using calibration constants:" << std::endl;
    for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
      try {
        fAmplCalibConst[p] = 1.2e-3 * fCalorimetryAlg.ElectronsFromADCPeak(1.0, p);
        std::cout << "   plane:" << p << " const:" << 1.0 / fAmplCalibConst[p] << std::endl;
      }
      catch (...) {
        fAmplCalibConst[p] = 1.0;
      }
    }
  }
  else {
    std::cout << "No plane-to-plane calibration." << std::endl;
    for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
      fAmplCalibConst[p] = 1.0;
    }
  }
  
  fDriftWindowInv = 1./fDriftWindow;
  
  std::string mode_str = "mean";
  std::cout  << "Downscale mode is: " << mode_str << std::endl;
  if (mode_str == "maxpool") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
    fDownscaleMode = kMax;
  }
  else if (mode_str == "maxmean") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMaxMean(dst, adc, tick0); };
    fDownscaleMode = kMaxMean;
  }
  else if (mode_str == "mean") {
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMean(dst, adc, tick0); };
    fDownscaleMode = kMean;
  }
  else {
    std::cout << "Downscale mode string not recognized, set to max pooling." << std::endl;
    //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
    fDownscaleMode = kMax;
  }
  
  
  //from https://internal.dunescience.org/doxygen/sp__fd__adcdump__job__example_8fcl_source.html
  /*
    fAdcMax    = 250.0;
    fAdcMin    = -5.0;
    Float_t OutMax  = 127.0; 
    Float_t OutMin = -128.0; 
  */ //anselmo
  fAdcMax    = 250.0;
  fAdcMin    = -5.0;
  Float_t OutMax  = 250.0; 
  Float_t OutMin =  -5; 
  
  
  fAdcMax    = 30.0;
  fAdcMin    = -10.0;
  OutMax  = 15.0; 
  OutMin =  -5.0; 
  
  
  
  fAdcOffset = OutMin;
  fAdcScale = (OutMax - fAdcOffset) / (fAdcMax - fAdcMin);
  fAdcZero   = fAdcOffset + fAdcScale * (0 - fAdcMin); // level of zero ADC after scaling;
  
  
  if (fAdcMax <= fAdcMin) {
    std::cout << "Misconfigured: AdcMax <= AdcMin" << std::endl;
  }
  if (fAdcScale == 0) {
    std::cout << "Misconfigured: OutMax == OutMin" << std::endl;
  }
  /*
    fBlurKernel = config.BlurKernel();
    fNoiseSigma = config.NoiseSigma();
    fCoherentSigma = config.CoherentSigma();
  */

  _tutils = new timeUtils(30);
#ifdef CompileTF
  fNNet->SetTimeUtils(_tutils);
#endif
}

CNNUtils::~CNNUtils(){

}

//*******************************************************
void CNNUtils::ComputeParticleCNN(AnaParticlePD& part){
//*******************************************************  


  
  std::cout << "CNN before: " << part.CNNscore[0] << " " << part.CNNscore[1] << " " << part.CNNscore[2] << " " << std::endl;

#ifdef CompileTF
  if (part.Hits[2].size()>0){
    for (size_t j=0;j<3;j++){
      part.CNNscore[j] = 0;
    }
    
    int nhits = 0;
    for(size_t ihit = 0; ihit < part.Hits[2].size(); ihit++){
      std::cout << " - " << ihit << ":" << part.Hits[2][ihit].CNN[0] << part.Hits[2][ihit].CNN[1] << part.Hits[2][ihit].CNN[2] << std::endl;
      for (size_t j=0;j<3;j++){
        part.CNNscore[j] += part.Hits[2][ihit].CNN[j];
      }
      nhits++;
    }
    
    if (nhits>=0){
      for (size_t j=0;j<3;j++){
        part.CNNscore[j] /= (Float_t)nhits;
      }
    }    
  }
#endif
  std::cout << "CNN after: " << part.CNNscore[0] << " " << part.CNNscore[1] << " " << part.CNNscore[2] << " " << std::endl;

}

//*******************************************************
void CNNUtils::produce(std::vector<AnaHitPD*>& hits, std::vector<AnaWireCNN>& wires){
//*******************************************************  

//  _tutils->printTime("time 0");

  _tutils->accumulateTime(0);

#ifndef CompileTF
  // sleep for 0.09 seconds to simulate the effect of TF when this is not enabled
  std::this_thread::sleep_for(std::chrono::milliseconds(9));
  _tutils->accumulateTime(1);
  //  return;  
#endif
  
  AnaHitPD* last_hit=NULL;
  std::vector< std::pair<unsigned int, float> > points;
  for (auto hit : hits){                  
    points.emplace_back(hit->WireID.Wire, hit->PeakTime);

    // 
    //    setWireDriftDataFromHit(*hit);
    last_hit=hit;
  }


  if (debug_levelA>=1){
    std::cout << " CNNUtils.cxx: produce. Processing hit: "  << std::endl;
    last_hit->Print();
  }
  
  size_t tpc = GetWireTPC(last_hit->WireID.Wire);

  setWireDriftData(wires,tpc, last_hit->WireID.Plane);
  
  _tutils->accumulateTime(1);
  
  //  _tutils->printTime("time 1");

  auto batch_out = predictIdVectors(points);
  if (debug_levelA>=1) std::cout << " CNNUtils.cxx: #points = "<<  points.size() << " batch_out.size() = "  << batch_out.size() << std::endl;

  if (points.size() != batch_out.size()){
      std::cout << "hits processing failed: " << points.size() << " " << batch_out.size()  << std::endl;
  }
  else{
    /*
    std::cout << last_hit->Channel << " " << last_hit->PeakTime << " ---> ";
    for (int jj=0;jj<batch_out[0].size();jj++)                                                                                                                                           
      std::cout << " " << batch_out[0][jj] << " ";                                                                                                                                       
    std::cout << std::endl;                                                                                                                                                              
    */
    if (last_hit)
      for (int jj=0;jj<3;jj++)
	last_hit->CNN[jj] = batch_out[0][jj];
    //    last_hit->Print();
    
  }

  _tutils->accumulateTime(2);
  //  _tutils->printTime("time 5");
}
  
//*******************************************************
std::vector<std::vector<float>> CNNUtils::predictIdVectors(std::vector<std::pair<unsigned int, float>> points){
//*******************************************************
  
  if (debug_levelA>=1) std::cout << wspacesA(2) << "predictIdVectors" << std::endl;   
  if (points.empty() || !fNNet) { return std::vector<std::vector<float>>(); }


  _tutils->accumulateTime(20);

  std::vector<std::vector<std::vector<float>>> inps(points.size(), std::vector<std::vector<float>>(fPatchSizeW, std::vector<float>(fPatchSizeD)));
  for (size_t i = 0; i < points.size(); ++i) {
    size_t wire = points[i].first;
    float drift = points[i].second;
    if (!bufferPatch(wire, drift, inps[i])) {std::cout << "Patch buffering failed" << std::endl;}    
  
  }
  _tutils->accumulateTime(21);

  std::vector<std::vector<float>> idvector = fNNet->Run(inps);

  _tutils->accumulateTime(22);

  return idvector;
}

//*******************************************************
bool CNNUtils::setWireDriftDataFromHit(const AnaHitPD& hit){
//*******************************************************

  if (debug_levelA>=1) std::cout << wspacesA(2) << "setWireDriftDataFromHit"<< std::endl;   

  size_t nwires = 480;//fGeometry->Nwires(plane, tpc, cryo);   anselmo
  size_t ndrifts = 6000;//det_prop.NumberTimeSamples();    anselmo

  _tutils->accumulateTime(10);

  // 1. 
  resizeView(nwires, ndrifts);

  _tutils->accumulateTime(11);
  size_t w_idx = hit.WireID.Wire;  

  std::vector<Float_t> adc(ndrifts,0);
  for (size_t i=0;i<hit.Signal.size();i++){
    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftDataFromHit. i, hit.Signal[i] = " << i << " " << hit.Signal[i] << " " << std::endl;   
    adc[i+hit.PeakTime-fPatchSizeD/2.] = hit.Signal[i];
  }
  
  if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftDataFromHit. w_idx, hit.Signal[i].size(), adc.size() = " << w_idx << " " << hit.Signal.size() << " " << adc.size()<< std::endl;   

  // 2. 
  auto wire_data = setWireData(adc, w_idx);
  _tutils->accumulateTime(12);
  //--------------------------
  size_t l=0;
  size_t non0_values=0;;
  size_t gt0_values=0;;
  for (auto v : wire_data) {
    if (v!=0){
      if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireDriftDataFromHit. Non 0 adc values in that channel. l, v = " << l << " " << v << std::endl;   
      non0_values++;
      if (v>0) gt0_values++;
    }
    l++;
  }

  if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftDataFromHit. wire_data.size() = " << wire_data.size() << std::endl;
  if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftDataFormHit. #non 0 adc values, >0 values = " << non0_values << " " << gt0_values << std::endl;   
  //--------------------------

  
  if (wire_data.empty()) {  // anselmo
    std::cout << "Wire data not set." << std::endl;
    return false; // also not critical, try to set other wires
  }
  fAlgView.fWireDriftData[w_idx] = wire_data;  // anselmo
  fAlgView.fWireChannels[w_idx] = hit.Channel;  
  
  //--------------------------  
  fAdcSumOverThr = 0;
  fAdcAreaOverThr = 0;

  for (auto v : adc) {
    if (v >= fAdcSumThr) {
      fAdcSumOverThr += v;
      fAdcAreaOverThr++;
    }
  }
  if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftDataFromHit. fAdcSumOverThr, fAdcSumAreaThr = " << fAdcSumOverThr << " " << fAdcAreaOverThr << std::endl;   


  if (debug_levelA>=10){ 
    for (auto v : wireData(w_idx)) {
      std::cout << wspacesA(8) << "setWireDriftDataFromHit. fAlgView.fWireDriftData[w_idx] = " << v << std::endl;   
    }
  }
  //--------------------------

      
  //applyBlur();
  //addWhiteNoise();
  //addCoherentNoise();
    
  return true;
}

//*******************************************************
bool CNNUtils::setWireDriftData(const std::vector<AnaWireCNN>& wires, Int_t hit_tpc, Int_t hit_plane){
//*******************************************************  


  // Set all wave forms for a given TPC and PLANE 

  
  if (debug_levelA>=1) std::cout << wspacesA(2) << "setWireDriftData"<< std::endl;   

  size_t nwires = 480;//fGeometry->Nwires(plane, tpc, cryo);   anselmo
  size_t ndrifts = 6000;//det_prop.NumberTimeSamples();    anselmo

  _tutils->accumulateTime(10);

  // 1. Give the appropriate size to the vectors and matrices used for storing the wave forms
  resizeView(nwires, ndrifts);

  _tutils->accumulateTime(11);

  if (debug_levelA>=2) std::cout << wspacesA(4) << "setWireDriftData. #wires = " << wires.size() << " hit tpc, plane = " << hit_tpc << " " << hit_plane << std::endl;   

  // LOOP OVER WIRES
  size_t iw=0;
  for (auto const& wire : wires) {  // loop over the wires

    size_t tpc   = GetWireTPC(wire.wire);
    size_t plane = GetWirePlane(wire.wire);

    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. wire = " << wire.wire << " --> tpc, plane = " << tpc << " " << plane << std::endl;   
    

    // Skip the were if it is not in the correct TPC and plane
    if ((Int_t)tpc != hit_tpc || (Int_t)plane != hit_plane){
      if (debug_levelA>=2) std::cout << wspacesA(6) << "setWireDriftData. i, wire, #adcs = " << iw << " " << wire.wire << " " << " " << wire.adcs.size()
                                     << " --> skip: tpc, plane = " << tpc << " " << plane << std::endl;         
      continue;
    }
    Int_t wireChannelNumber=0;
    //    auto wireChannelNumber = wire.Channel();  // TODO. Get the cannel number for that wire

    // Need a wire index between 0 and 479
    size_t w_idx = GetWireIndex(wire.wire);

    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. w_idx = " << w_idx << std::endl;   
    
    // Fill the adc vector 
    std::vector<Float_t> adc(ndrifts,0);
    for (size_t i=0;i<wire.adcs.size();i++){
      if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireDriftData. i, adcs[i] = " << i << " " << wire.adcs[i] << " " << std::endl;   
      adc[wire.time+i] = wire.adcs[i];
    }
    
    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. w_idx, adcs[i].size(), adc.size() = " << w_idx << " " << wire.adcs.size() << " " << adc.size()<< std::endl;   


    
    if (debug_levelA>=2) std::cout << wspacesA(6) << "setWireDriftData. i, wire, w_idx, #adcs = " << iw << " " << wire.wire << " " << w_idx << " " << wire.adcs.size() << std::endl;   
    iw++;
    
    // 2. 
    auto wire_data = setWireData(adc, w_idx);
    _tutils->accumulateTime(12);
    //--------------------------
    size_t l=0;
    size_t non0_values=0;;
    size_t gt0_values=0;;
    for (auto v : wire_data) {
      if (v!=0){
        if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireDriftData. Non 0 adc values in that channel. l, v = " << l << " " << v << std::endl;   
        non0_values++;
        if (v>0) gt0_values++;
      }
      l++;
    }
    
    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. wire_data.size() = " << wire_data.size() << std::endl;
    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. #non 0 adc values, >0 values = " << non0_values << " " << gt0_values << std::endl;   
    //--------------------------
    
    
    if (wire_data.empty()) {  // anselmo
      std::cout << "Wire data not set." << std::endl;
      continue; // also not critical, try to set other wires
    }
    fAlgView.fWireDriftData[w_idx] = wire_data;  
    fAlgView.fWireChannels[w_idx]  = wireChannelNumber;    // NOT REALLY USED
    
    //--------------------------  
    fAdcSumOverThr = 0;
    fAdcAreaOverThr = 0;
    
    for (auto v : adc) {
      if (v >= fAdcSumThr) {
        fAdcSumOverThr += v;
        fAdcAreaOverThr++;
      }
    }
    if (debug_levelA>=3) std::cout << wspacesA(6) << "setWireDriftData. fAdcSumOverThr, fAdcSumAreaThr = " << fAdcSumOverThr << " " << fAdcAreaOverThr << std::endl;   


    if (debug_levelA>=10){ 
      for (auto v : wireData(w_idx)) {
        std::cout << wspacesA(8) << "setWireDriftData. fAlgView.fWireDriftData[w_idx] = " << v << std::endl;   
      }
    }    
  }

  //--------------------------

      
  //applyBlur();
  //addWhiteNoise();
  //addCoherentNoise();
    
  return true;
}

//*******************************************************
void CNNUtils::resizeView(size_t wires,size_t drifts){
//*******************************************************

  /*
    Resize the vectors and matrices 

     - fWireChannels  (#wires)
     - fWireDriftData  (#wires x #ChachedDrifts)
   */

  
  _tutils->accumulateTime(15);

  if (debug_levelA>=2) std::cout << wspacesA(4) << "resizeView. total drifts, total wires = " << drifts << " " << wires << std::endl;   
  DataProviderAlgView& result = fAlgView;
  result.fNWires = wires;
  result.fNDrifts = drifts;
  result.fNScaledDrifts = drifts / fDriftWindow;
  result.fNCachedDrifts = fDownscaleFullView ? result.fNScaledDrifts : drifts;

  _tutils->accumulateTime(16);

  //  result.fWireChannels.resize(wires, raw::InvalidChannelID);
  result.fWireChannels.resize(wires, 4294967295); 
  result.fWireDriftData.resize(wires, std::vector<float>(result.fNCachedDrifts, fAdcZero));   

  if (debug_levelA>=2) std::cout << wspacesA(4) << "resizeView.  final fWireDriftData (NxM) = " << wires << " x " << result.fNCachedDrifts << std::endl;   
  
  _tutils->accumulateTime(17);

  if (fCalibrateLifetime) {
    result.fLifetimeCorrFactors.resize(drifts);
    for (size_t t = 0; t < drifts; ++t) {
      result.fLifetimeCorrFactors[t] = fCalorimetryAlg.LifetimeCorrection(t);
    }
  }
  /*
  else {
    for (size_t t = 0; t < drifts; ++t) {
      result.fLifetimeCorrFactors[t] = 1.0;
    }
  }
  */

  _tutils->accumulateTime(18);

  //  return result;
}


//*******************************************************
std::vector<float> CNNUtils::setWireData(std::vector<float> const& adc, size_t wireIdx) const{
//*******************************************************

  /* RETURNS THE ADC VECTOR FOR A GIVEN WIRE, BUT DOES NOT SET IT. SEVERAL OPTIONS
     1. wireIndx > fAlgView.fWireDriftData.size() --> returns empty vector
     2. adc.empty()                               --> returns empty vector; 
     3. fDownscaleFullView = true                 --> returns a downscaled version of the input adc vector
     4. fDownscaleFullView = false: 
        4.1. adc.size() <= wData.size()                --> returns input adc vector
        4.2. adc.size() > wData.size()                 --> returns first wData.size() from input adc vector  
   */

  if (debug_levelA>=3) std::cout << wspacesA(4) << "setWireData. wireIdx, fAlgView.fWireDriftData.size() = " << wireIdx << " " << fAlgView.fWireDriftData.size() << std::endl;   

  // Nothing to do if the wire index is above the size of the vector
  if (wireIdx >= fAlgView.fWireDriftData.size()){
    if (debug_levelA>=3) std::cout << wspacesA(4) << "setWireData. wireIdx >= fAlgView.fWireDriftData.size() --> THERE IS A PROBLEM !!!"<< std::endl;
    return std::vector<float>();//std::nullopt;  
  }

  if (debug_levelA>=3) std::cout << wspacesA(4) << "setWireData --> OK !!!"<< std::endl;
  auto& wData = fAlgView.fWireDriftData[wireIdx];
  
  if (fDownscaleFullView) {
    if (!adc.empty()) { return downscale(wData.size(), adc, 0); }
    else  {return std::vector<float>();}//std::nullopt;}
  }
  else {
    if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireData. adc.size(). fWireDriftData[wireIdx].size() = " << adc.size() << " " << wData.size()<< std::endl;             
    if (adc.empty()) { std::vector<float>();}//std::nullopt; }
    else if (adc.size() <= wData.size()){
      // Returns the entire input vector
      if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireData. 1. adc.size() " << adc.size() << std::endl;
      return adc;
    }
    else {
      // Returns only a fraction of the input vector
      if (debug_levelA>=4) std::cout << wspacesA(8) << "setWireData. 2. adc.size() " << adc.size() << std::endl;                
      return std::vector<float>(adc.begin(), adc.begin() + wData.size());
    }
  }
  //  return std::make_optional(wData);
  return wData;
}



//*******************************************************
bool CNNUtils::bufferPatch(size_t wire, float drift, std::vector<std::vector<float>>& patch) const {
//*******************************************************
  
  if (debug_levelA>=2) std::cout << wspacesA(4) << "bufferPatch: wire, drift, patch.size() = " << wire << " " << drift << " " << patch.size() << std::endl;   
  if (fDownscaleFullView) {
    size_t sd = (size_t)(drift / fDriftWindow);
    if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == sd))
      return true; // still within the current position
    
    fCurrentWireIdx = wire;
    fCurrentScaledDrift = sd;
    
    return patchFromDownsampledView(wire, drift, fPatchSizeW, fPatchSizeD, patch);
  }
  else {
    if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == drift))
      return true; // still within the current position
    
    fCurrentWireIdx = wire;
    fCurrentScaledDrift = drift;
    
    return patchFromOriginalView(wire, drift, fPatchSizeW, fPatchSizeD, patch);
  }
}

//*******************************************************
bool CNNUtils::patchFromOriginalView(size_t wire,float drift,size_t size_w,size_t size_d,std::vector<std::vector<float>>& patch) const{
//*******************************************************

  /*


   */
  
  int dsize = fDriftWindow * size_d;
  int halfSizeW = size_w / 2;
  int halfSizeD = dsize / 2;
  
  int w0 = wire - halfSizeW;
  int w1 = wire + halfSizeW;
  
  int d0 = int(drift) - halfSizeD;
  int d1 = int(drift) + halfSizeD;
  
  if (d0 < 0) d0 = 0;
  
  if (debug_levelA>=2) std::cout << wspacesA(6) << "patchFromOriginalView. wire, drift, w0-w1, d0-d1 = " << wire << " " << drift << " " << w0 << "-" << w1 << " " << d0 << "-" << d1 << std::endl;   
  
  std::vector<float> adc_down(dsize);
  int wsize = fAlgView.fWireDriftData.size();
  for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch) {
    if (debug_levelA>=5) std::cout << wspacesA(10) << "patchFromOriginalView. w, wpatch = " << w << " " << wpatch << std::endl;   
    size_t count_gt0=0;  // anselmo
    if ((w >= 0) && (w < wsize)) {
      auto& src = fAlgView.fWireDriftData[w];
      if (debug_levelA>=6) std::cout << wspacesA(12) << "patchFromOriginalView. fWireDriftData[w].size() = " << fAlgView.fWireDriftData[w].size() << std::endl;   
      int src_size = src.size();
      for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch) {
        if ((d >= 0) && (d < src_size)) { adc_down[dpatch] = src[d]; }
        else {
          adc_down[dpatch] = fAdcZero;
        }
        if (debug_levelA>=7) std::cout << wspacesA(14) << "patchFromOriginalView. d,dpatch,src[d]  = " << d << " " << dpatch << " " << src[d] << std::endl;   
        if (src[d]>0) count_gt0++;  // anselmo
      }
    }
    else {
      std::fill(adc_down.begin(), adc_down.end(), fAdcZero);
    }

    if (debug_levelA>=2 && ( count_gt0>0 || debug_levelA>=4)){
      std::cout << wspacesA(8) << "patchFromOriginalView. After downscale: w, wpatch, adc_down.size(), patch[wpatch].size(), #>0 = "
                << w << " "
                << wpatch << " " 
                << adc_down.size() << " "
                << patch[wpatch].size()
                << " " << count_gt0 << std::endl;   
    }
      
    patch[wpatch] = downscale(patch[wpatch].size(), adc_down, d0, count_gt0);

  }

  if (debug_levelA>=5) std::cout << wspacesA(10) << "patchFromOriginalView. patch.size() = " << patch.size() << std::endl;   
  return true;
}


//*******************************************************
std::vector<float> CNNUtils::downscale(std::size_t dst_size, std::vector<float> const& adc, size_t tick0, size_t adc_gt0) const{    
//*******************************************************

  if (debug_levelA>=2 && (adc_gt0 || debug_levelA>=4)) std::cout << wspacesA(10) << "downscale. mode, dst_size,  " << fDownscaleMode << " " << dst_size << std::endl;   

  switch (fDownscaleMode) {
  case kMean:    return downscaleMean(dst_size, adc, tick0, adc_gt0);
  case kMaxMean: return downscaleMaxMean(dst_size, adc, tick0);
  case kMax:     return downscaleMax(dst_size, adc, tick0);
  }
  std::cout << "Downscale mode not supported." << std::endl;
  return std::vector<float>();
}

//*******************************************************
std::vector<float> CNNUtils::downscaleMean(std::size_t dst_size,std::vector<float> const& adc,size_t tick0, size_t adc_gt0) const{
//*******************************************************
  
  if (debug_levelA>=2 && (adc_gt0>0 || debug_levelA>=4)) std::cout << wspacesA(12) << "downscaleMean. adc.size(), tick0 = " << adc.size() << " " << tick0 << std::endl;   
  size_t count=0; //anselmo

  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }

  // loop over wires in groups of fDriftWindow batches 
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;
    if (debug_levelA>=8) std::cout << wspacesA(16) << "downscaleMean, i = " << i << std::endl;   

    // sum all adc values in a batch applying Lifetime correction
    float sum_adc = 0;
    for (size_t k = k0; k < k1; ++k) {
      if (fCalibrateLifetime){
	if (k + tick0 < fAlgView.fLifetimeCorrFactors.size())
	  sum_adc += adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      }
      else
	  sum_adc += adc[k];

      if (debug_levelA>=8) std::cout << wspacesA(16) << "downscaleMean. adc["<< k << "] = " << adc[k] << std::endl;   
    }
    result[i] = sum_adc * fDriftWindowInv;
    
    if (result[i]>0){
      if (debug_levelA>=8) std::cout << wspacesA(16) << "downscaleMean. >0 result["<< i << "] = " << result[i] << std::endl;   
      count++;
    }
  }
  if (debug_levelA>=7) std::cout << wspacesA(14) << "downscaleMean. result.size(), >0 results = " << result.size() << " " << count << std::endl;   
      scaleAdcSamples(result, adc_gt0);

  count=0; // anselmo
  if (debug_levelA>=7) std::cout << wspacesA(14) << "downscaleMean. results after downscaling. result.size() = " << result.size() << std::endl;   
  for (size_t i=0;i<result.size();i++){
    if (debug_levelA>=8) std::cout << wspacesA(16) << "downscaleMean. result["<<i<<"] = " << result[i] << std::endl;   
    if (result[i]>0){ count++;
      if (debug_levelA>=3 && (adc_gt0>0 || debug_levelA>=4)) std::cout << wspacesA(12) << "downscaleMean. result["<<i<<"] = " << result[i] << std::endl;      
    }
  }
  if (debug_levelA>=2 && (adc_gt0>0 || debug_levelA>=4)) std::cout << wspacesA(12) << "downscaleMean. results after downscaling ---->  > 0 results = " << count << std::endl;   
    
  return result;
}

//*******************************************************
    void CNNUtils::scaleAdcSamples(std::vector<float>& values, size_t adc_gt0) const{
//*******************************************************
  
      if (debug_levelA>=2 && (adc_gt0>0 || debug_levelA>=4)) std::cout << wspacesA(14) << "scaleAdcSamples. fPlane = " << fPlane << std::endl;   
  float calib = fAmplCalibConst[fPlane];
  auto* data = values.data();
  size_t k = 0, size4 = values.size() >> 2, size = values.size();
  if (debug_levelA>=8) std::cout << wspacesA(16) << "scaleAdcSamples. size4, size, calib = " << size4 << " " << size << " " << calib << std::endl;   
  for (size_t i = 0; i < size4; ++i) // vectorize if you can
    {
      
      data[k] *= calib; // prescale by plane-to-plane calibration factors
      data[k + 1] *= calib;
      data[k + 2] *= calib;
      data[k + 3] *= calib;

      if (data[k] < fAdcMin) { data[k] = fAdcMin; } // saturate min
      if (data[k + 1] < fAdcMin) { data[k + 1] = fAdcMin; }
      if (data[k + 2] < fAdcMin) { data[k + 2] = fAdcMin; }
      if (data[k + 3] < fAdcMin) { data[k + 3] = fAdcMin; }
      
      if (data[k] > fAdcMax) { data[k] = fAdcMax; } // saturate max
      if (data[k + 1] > fAdcMax) { data[k + 1] = fAdcMax; }
      if (data[k + 2] > fAdcMax) { data[k + 2] = fAdcMax; }
      if (data[k + 3] > fAdcMax) { data[k + 3] = fAdcMax; }

      data[k] = fAdcOffset +
	fAdcScale *
	(data[k] - fAdcMin); // shift and scale to the output range, shift to the output min
      data[k + 1] = fAdcOffset + fAdcScale * (data[k + 1] - fAdcMin);
      data[k + 2] = fAdcOffset + fAdcScale * (data[k + 2] - fAdcMin);
      data[k + 3] = fAdcOffset + fAdcScale * (data[k + 3] - fAdcMin);

      k += 4;
      if (debug_levelA>=9) std::cout << wspacesA(18) << "scaleAdcSamples. i, k, data[k], values[k] = " << i << " " << k << " " << data[k] << " " << values[k] << std::endl;   
    }
  while (k < size) {
    if (debug_levelA>=9) std::cout << wspacesA(18) << "scaleAdcSamples. k, data[k], values[k] = " << k << " " << data[k] << " " << values[k] << std::endl;   
    data[k] = scaleAdcSample(data[k]);
    ++k;
  }// do the tail
}


//*******************************************************
float CNNUtils::scaleAdcSample(float val) const{
//*******************************************************  

  if (debug_levelA>=2) std::cout << wspacesA(16) << "scaleAdcSample" << std::endl;   
  if (debug_levelA>=7) std::cout << wspacesA(16) << "scaleAdcSample. val, Plane, fAmplCalibConst[fPlane] = " << val << " " << fPlane << " " << fAmplCalibConst[fPlane] << std::endl;   
  
  val *= fAmplCalibConst[fPlane]; // prescale by plane-to-plane calibration factors
  
  if (val < fAdcMin) { val = fAdcMin; } // saturate
  else if (val > fAdcMax) {
    val = fAdcMax;
  }
  
  return fAdcOffset + fAdcScale * (val - fAdcMin); // shift and scale to the output range, shift to the output min
}

// ------------------------------------------------------
// -------------------ModelInterface---------------------
// ------------------------------------------------------

//*******************************************************
std::vector<std::vector<float>> nnet2::ModelInterface::Run(std::vector<std::vector<std::vector<float>>> const& inps, int samples) {
//*******************************************************
  
  if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty())
    return std::vector<std::vector<float>>();

  if ((samples == -1) || (samples > (int)inps.size())) { samples = inps.size(); }

  std::vector<std::vector<float>> results;
  for (int i = 0; i < samples; ++i) {
    results.push_back(Run(inps[i]));
  }
  return results;
}

//*******************************************************
std::string nnet2::ModelInterface::findFile(const char* fileName) const{
//*******************************************************

  std::string fname_out;
  std::string path = "";
  /*
  cet::search_path sp("FW_SEARCH_PATH");
  if (!sp.find_file(fileName, fname_out)) {
    struct stat buffer;
    if (stat(fileName, &buffer) == 0) { fname_out = fileName; }
    else {
      throw art::Exception(art::errors::NotFound) << "Could not find the model file " << fileName;
    }
  }
  */
  fname_out =  path + fileName;
  std::cout << fname_out << std::endl;
  return fname_out;
}


// ------------------------------------------------------
// -----------------TfModelInterface---------------------
// ------------------------------------------------------
//*******************************************************
nnet2::TfModelInterface::TfModelInterface(const char* modelFileName){
//*******************************************************  

#ifdef  CompileTF
  g = tf::Graph::create(nnet2::ModelInterface::findFile(modelFileName).c_str(),
                        {"cnn_output", "_netout"});
  //if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
  if (!g) {std::cout << "TF model failed." << std::endl; exit(1);}
#endif
  std::cout << "TF model loaded." << std::endl;
  
}
// ------------------------------------------------------

//*******************************************************
std::vector<std::vector<float>> nnet2::TfModelInterface::Run(std::vector<std::vector<std::vector<float>>> const& inps, int samples){
//*******************************************************
  
  if (debug_levelA>=2) std::cout << wspacesA(4) << "Run. inps.size() = " << inps.size() << std::endl;   
  if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty())
    return std::vector<std::vector<float>>();

  if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

  long long int rows = inps.front().size(), cols = inps.front().front().size();

  if (debug_levelA>=3) std::cout << wspacesA(6) << "Run. rows, cols = " << rows << " " << cols << std::endl;   
  
  if (debug_levelA>=3) std::cout << wspacesA(6) << "Run. READY TO CALL TF " << std::endl;   

  size_t count_gt0=0;

#ifdef  CompileTF
  tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({samples, rows, cols, 1}));   // anselmo
  auto input_map = _x.tensor<float, 4>();
#endif
  for (long long int s = 0; s < samples; ++s) {
    const auto& sample = inps[s];
    for (long long int r = 0; r < rows; ++r) {
      const auto& row = sample[r];
      for (long long int c = 0; c < cols; ++c) {
#ifdef  CompileTF
        input_map(s, r, c, 0) = row[c];  // anselmo
#endif
        if (row[c]>0) count_gt0++;
        if (debug_levelA>=5) std::cout << wspacesA(8) << "Run. s, r, c, row[c] = " << s << " " << r << " " << c << " " << row[c] << std::endl;
        else if (debug_levelA>=4 && row[c]>0 ) std::cout << wspacesA(6) << "Run. s, r, c, row[c] = " << s << " " << r << " " << c << " " << row[c] << std::endl;   
      }
    }
  }

  if (debug_levelA>=3) std::cout << wspacesA(6) << "Run. #tensor elements > 0 = " << count_gt0 << std::endl;   

#ifdef  CompileTF 
  return g->run(_x);
#else
  return std::vector<std::vector<float>>();
#endif
}
// ------------------------------------------------------

//*******************************************************
std::vector<float> nnet2::TfModelInterface::Run(std::vector<std::vector<float>> const& inp2d){
//*******************************************************
  
  long long int rows = inp2d.size(), cols = inp2d.front().size();

  if (debug_levelA>=2) std::cout << wspacesA(4) << "Run. rows, cols = " << rows << " " << cols << std::endl;   
  /*
  tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, rows, cols, 1}));
  auto input_map = _x.tensor<float, 4>();
  for (long long int r = 0; r < rows; ++r) {
    const auto& row = inp2d[r];
    for (long long int c = 0; c < cols; ++c) {
      input_map(0, r, c, 0) = row[c];
    }
  }

  auto out = g->run(_x);
  if (!out.empty())
    return out.front();
  else
    return std::vector<float>();
  */
  return std::vector<float>();
}



//*******************************************************
size_t CNNUtils::GetWireTPC(Int_t wire){
//*******************************************************

  Int_t tpcs[6]={1,5,9,2,6,10};

  for (size_t i = 0;i<6;i++){
    Int_t tpc = tpcs[i];
    for (size_t plane = 0;plane<3;plane++){
      if (wire>wire_mapping[tpc][plane].first && wire<=wire_mapping[tpc][plane].second)
        return tpc;       
    }
  }
  return 0;
}

//*******************************************************
size_t CNNUtils::GetWirePlane(Int_t wire){
//*******************************************************


  Int_t tpcs[6]={1,5,9,2,6,10};

  for (size_t i = 0;i<6;i++){
    Int_t tpc = tpcs[i];
    for (size_t plane = 0;plane<3;plane++){
      if (wire>wire_mapping[tpc][plane].first && wire<=wire_mapping[tpc][plane].second)
        return plane;       
    }
  }
  return 0;

}


//*******************************************************
size_t CNNUtils::GetWireIndex(Int_t wire){
//*******************************************************

  Int_t tpcs[6]={1,5,9,2,6,10};

  for (size_t i = 0;i<6;i++){
    Int_t tpc = tpcs[i];
    for (size_t plane = 0;plane<3;plane++){
      if (wire>wire_mapping[tpc][plane].first && wire<=wire_mapping[tpc][plane].second)
        return wire-wire_mapping[tpc][plane].first;       
    }
  }
  return 0;
  
}

