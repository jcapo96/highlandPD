#include "CNNUtils.hxx"

//*******************************************************
void CNNUtils::produce(std::vector<AnaHitPD>& hits){
//*******************************************************  

  std::vector< std::pair<unsigned int, float> > points;
  for (auto const & hit : hits){                  
    points.emplace_back(hit.WireID.Wire, hit.PeakTime);
    fPointIdAlg.setWireDriftDataFromHit(hit);
    std::cout << hit.Channel << " " << hit.PeakTime << " ---> ";
  }
  auto batch_out = fPointIdAlg.predictIdVectors(points);
  if (debug_level>=1) std::cout << " CNNUtils.cxx: #points = "<<  points.size() << " batch_out.size() = "  << batch_out.size() << std::endl;
  if (points.size() != batch_out.size()){
      //                        throw cet::exception("CNNUtils") << "hits processing failed" << std::endl;
      std::cout << "hits processing failed: " << points.size() << " " << batch_out.size()  << std::endl;
  }

																						       
  for (int jj=0;jj<batch_out[0].size();jj++)                                                                                                                                           
    std::cout << " " << batch_out[0][jj] << " ";                                                                                                                                       
  std::cout << std::endl;                                                                                                                                                              

}
  
//*******************************************************
void CNNUtils::FillHits(const AnaEventPD& evt, std::vector<AnaHitPD*>& hitPtrList){
//*******************************************************

  
  if (debug_level>=0) std::cout << "CNNUtils.cxx: #parts "<< evt.nParticles << std::endl;  
  for (Int_t i=0;i<evt.nParticles/2+1;i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(evt.Particles[i]);

    std::cout << "ANSELMO: particle " << i << " --> #hits = " << part->Hits[2].size() << std::endl;       
    
    if (debug_level>=0){
      std::cout << "CNNUtils.cxx: hits list for this particle: " << std::endl;       
      for (UInt_t j=0;j<part->Hits[2].size();j++){

        AnaHitPD& hit = part->Hits[2][j];
        
        std::cout << " - " << j << ": wire, time, ampl, integral =  "
                  << hit.WireID.Wire << " "
                  << hit.PeakTime << " "
                  << hit.PeakAmplitude << " "
                  << hit.Integral << " "
                  << hit.Channel << " "
                  << hit.StartTick << "-"
                  << hit.EndTick << " "
          //                  << evt.ADC[hit.Channel][hit.StartTick] << " "
                  << hit.Signal.size() << " " 
                  << std::endl;       
      }
    }
    
    for (UInt_t j=0;j<part->Hits[2].size();j++){

      AnaHitPD& hit = part->Hits[2][j];

      // Fill hits with some "meaningful"  adc values
      /*
      for (size_t k=hit.fPeakTime-5;k<hit.fPeakTime+5;k++){
        hit.fSignal.push_back(hit.fPeakAmplitude*10);
      }
      */
      fPointIdAlg.setWireDriftDataFromHit(hit);
      //      hitListHandle.push_back(&(hit));
      hitPtrList.push_back(&(hit));
    }
  }
}

//*******************************************************
bool img::DataProviderAlg::setWireDriftDataFromHit(const AnaHitPD& hit){
//*******************************************************

  if (debug_level>=1) std::cout << wspaces(2) << "DataProviderAlg::setWireDriftDataFromHit"<< std::endl;   

  static bool first=true;

  
  fAdcSumOverThr = 0;
  fAdcAreaOverThr = 0;
  
  size_t nwires = 480;//fGeometry->Nwires(plane, tpc, cryo);   anselmo
  size_t ndrifts = 6000;//det_prop.NumberTimeSamples();    anselmo

  fAlgView = resizeView(clock_data, det_prop, nwires, ndrifts);
    
  size_t w_idx = hit.WireID.Wire;  

  size_t fPatchSizeD = 48;
  std::vector<Float_t> adc(ndrifts,0);
  for (size_t i=0;i<hit.Signal.size();i++)
    adc[i+hit.PeakTime-fPatchSizeD/2.] = hit.Signal[i];
  
  
  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftDataFromHit. Found wireID in required TPC,cryo and plane. w_idx, hit.Signal[i].size(), adc.size() = " << w_idx << " " << hit.Signal.size() << " " << adc.size()<< std::endl;   
  auto wire_data = setWireData(adc, w_idx);
  
  size_t l=0;
  size_t non0_values=0;;
  size_t gt0_values=0;;
  for (auto v : wire_data) {
    if (v!=0){
      if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireDriftDataFromHit. Non 0 adc values in that channel. l, v = " << l << " " << v << std::endl;   
      non0_values++;
      if (v>0) gt0_values++;
    }
    l++;
  }

  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftDataFromHit. wire_data.size() = " << wire_data.size() << std::endl;
  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftDataFormHit. #non 0 adc values, >0 values = " << non0_values << " " << gt0_values << std::endl;   
  //    if (!wire_data) {
  if (wire_data.empty()) {  // anselmo
    //mf::LogWarning("DataProviderAlg") << "Wire data not set.";
    std::cout << "Wire data not set." << std::endl;
    return false; // also not critical, try to set other wires
  }
  //    fAlgView.fWireDriftData[w_idx] = *wire_data;
  fAlgView.fWireDriftData[w_idx] = wire_data;  // anselmo
  for (auto v : adc) {
    if (v >= fAdcSumThr) {
      fAdcSumOverThr += v;
      fAdcAreaOverThr++;
    }
  }
  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftDataFromHit. fAdcSumOverThr, fAdcSumAreaThr = " << fAdcSumOverThr << " " << fAdcAreaOverThr << std::endl;   
  
  fAlgView.fWireChannels[w_idx] = hit.Channel;

  if (debug_level>=10) 
    for (auto v : wireData(w_idx)) {
      std::cout << wspaces(8) << "DataProviderAlg::setWireDriftDataFromHit. fAlgView.fWireDriftData[w_idx] = " << v << std::endl;   
    }

  
  //applyBlur();
  //addWhiteNoise();
  //addCoherentNoise();
    
  return true;
}

img::DataProviderAlgView
img::DataProviderAlg::resizeView(size_t wires,size_t drifts){

  if (debug_level>=2) std::cout << wspaces(4) << "DataProviderAlg::resizeView. total drifts, total wires = " << drifts << " " << wires << std::endl;   
  img::DataProviderAlgView result;
  result.fNWires = wires;
  result.fNDrifts = drifts;
  result.fNScaledDrifts = drifts / fDriftWindow;
  result.fNCachedDrifts = fDownscaleFullView ? result.fNScaledDrifts : drifts;
  
  //  result.fWireChannels.resize(wires, raw::InvalidChannelID);
  result.fWireChannels.resize(wires, 4294967295);


  
  result.fWireDriftData.resize(wires, std::vector<float>(result.fNCachedDrifts, fAdcZero)); 
  
  result.fLifetimeCorrFactors.resize(drifts);
  if (fCalibrateLifetime) {
    for (size_t t = 0; t < drifts; ++t) {
      //      result.fLifetimeCorrFactors[t] = fCalorimetryAlg.LifetimeCorrection(clock_data, det_prop, t);
      //     std::cout << "DataProviderAlg::LifeTimeCorrection" << std::endl;   
      result.fLifetimeCorrFactors[t] = fCalorimetryAlg.LifetimeCorrection(t);
    }
  }
  else {
    for (size_t t = 0; t < drifts; ++t) {
      result.fLifetimeCorrFactors[t] = 1.0;
    }
  }
  return result;
}
