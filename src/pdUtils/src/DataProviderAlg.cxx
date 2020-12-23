
// Class:       PointIdAlg
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
     
//#include "larreco/RecoAlg/ImagePatternAlgs/DataProvider/DataProviderAlg.h"
#include "DataProviderAlg.hxx"


//#include "art/Framework/Services/Registry/ServiceHandle.h"
     
//#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
   
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
//#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
//#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

//#include "messagefacility/MessageLogger/MessageLogger.h"
namespace detinfo {
  class DetectorProperties;
}
namespace geo {
  class GeometryCore;
}
    
//#include "CLHEP/Random/RandGauss.h"
    
/*img::DataProviderAlg::DataProviderAlg(const Config& config)
  : fAlgView{}
  , fDownscaleMode(img::DataProviderAlg::kMax)
  , fDriftWindow(10)
  , fCalorimetryAlg(config.CalorimetryAlg())
  , fGeometry(art::ServiceHandle<geo::Geometry const>().get())
  , fAdcSumOverThr(0)
  , fAdcSumThr(10)
  , // set fixed threshold of 10 ADC counts for counting the sum
  fAdcAreaOverThr(0)
  , fNoiseSigma(0)
  , fCoherentSigma(0)
  {
    fCalibrateLifetime = config.CalibrateLifetime();
    fCalibrateAmpl = config.CalibrateAmpl();
    
    fAmplCalibConst.resize(fGeometry->MaxPlanes());
    if (fCalibrateAmpl) {
      mf::LogInfo("DataProviderAlg") << "Using calibration constants:";
      for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
	try {
          fAmplCalibConst[p] = 1.2e-3 * fCalorimetryAlg.ElectronsFromADCPeak(1.0, p);
          mf::LogInfo("DataProviderAlg") << "   plane:" << p << " const:" << 1.0 / fAmplCalibConst[p];
        }
	catch (...) {
          fAmplCalibConst[p] = 1.0;
        }
      }
    }
    else {
      mf::LogInfo("DataProviderAlg") << "No plane-to-plane calibration.";
      for (size_t p = 0; p < fAmplCalibConst.size(); ++p) {
	fAmplCalibConst[p] = 1.0;
      }
    }
    
    fDriftWindow = config.DriftWindow();
    fDownscaleFullView = config.DownscaleFullView();
    fDriftWindowInv = 1.0 / fDriftWindow;
    
    std::string mode_str = config.DownscaleFn();
    mf::LogVerbatim("DataProviderAlg") << "Downscale mode is: " << mode_str;
    if (mode_str == "maxpool") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMax;
    }
    else if (mode_str == "maxmean") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMaxMean(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMaxMean;
    }
    else if (mode_str == "mean") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMean(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMean;
    }
    else {
      mf::LogError("DataProviderAlg") << "Downscale mode string not recognized, set to max pooling.";
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMax;
    }
    
    fAdcMax = config.AdcMax();
    fAdcMin = config.AdcMin();
    fAdcOffset = config.OutMin();
    fAdcScale = (config.OutMax() - fAdcOffset) / (fAdcMax - fAdcMin);
    fAdcZero = fAdcOffset + fAdcScale * (0 - fAdcMin); // level of zero ADC after scaling
    
    if (fAdcMax <= fAdcMin) {
      throw cet::exception("img::DataProviderAlg") << "Misconfigured: AdcMax <= AdcMin" << std::endl;
    }x
    if (fAdcScale == 0) {
      throw cet::exception("img::DataProviderAlg") << "Misconfigured: OutMax == OutMin" << std::endl;
    }
    
    fBlurKernel = config.BlurKernel();
    fNoiseSigma = config.NoiseSigma();
    fCoherentSigma = config.CoherentSigma();
  }
// ------------------------------------------------------
*/



std::string wspaces(size_t n){
  std::string ws="";
  for (size_t i=0;i<n;i++)
    ws +=" ";
  return ws;
}

img::DataProviderAlg::DataProviderAlg()
  : fAlgView{}
  , fDownscaleMode(img::DataProviderAlg::kMax)
  , fDriftWindow(10)
    //  , fCalorimetryAlg(config.CalorimetryAlg())
    //  , fGeometry(art::ServiceHandle<geo::Geometry const>().get())
  , fAdcSumOverThr(0)
  , fAdcSumThr(10)
  , // set fixed threshold of 10 ADC counts for counting the sum
  fAdcAreaOverThr(0)
    //  , fNoiseSigma(0)
    //  , fCoherentSigma(0)
  {
    fCalibrateLifetime = true;
    fCalibrateAmpl = true;
    
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

    fDriftWindow = 6;
    fDownscaleFullView = false;
    fDriftWindowInv = 1./fDriftWindow;
    
    std::string mode_str = "mean";
    std::cout  << "Downscale mode is: " << mode_str << std::endl;
    if (mode_str == "maxpool") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMax;
    }
    else if (mode_str == "maxmean") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMaxMean(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMaxMean;
    }
    else if (mode_str == "mean") {
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMean(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMean;
    }
    else {
      std::cout << "Downscale mode string not recognized, set to max pooling." << std::endl;
      //fnDownscale = [this](std::vector<float> & dst, std::vector<float> const & adc, size_t tick0) { downscaleMax(dst, adc, tick0); };
      fDownscaleMode = img::DataProviderAlg::kMax;
    }


    //from https://internal.dunescience.org/doxygen/sp__fd__adcdump__job__example_8fcl_source.html
    fAdcMax    = 250.0;
    fAdcMin    = -5.0;
    Float_t OutMax  = 127.0; 
    Float_t OutMin = -128.0; 
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
  }


// ------------------------------------------------------
img::DataProviderAlg::~DataProviderAlg() = default;
// ------------------------------------------------------

img::DataProviderAlgView
img::DataProviderAlg::resizeView(detinfo::DetectorClocksData const& clock_data,
				 detinfo::DetectorPropertiesData const& det_prop,
				 size_t wires,
				 size_t drifts)
{

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
// ------------------------------------------------------

float
img::DataProviderAlg::poolMax(int wire, int drift, size_t r) const
{
  std::cout << "DataProviderAlg::poolMax" << std::endl;   
  size_t rw = r, rd = r;
  if (!fDownscaleFullView) { rd *= fDriftWindow; }
  
  size_t didx = getDriftIndex(drift);
  int d0 = didx - rd;
  if (d0 < 0) { d0 = 0; }
  int d1 = didx + rd;
  if (d1 >= (int)fAlgView.fNCachedDrifts) { d1 = fAlgView.fNCachedDrifts - 1; }
  
  int w0 = wire - rw;
  if (w0 < 0) { w0 = 0; }
  int w1 = wire + rw;
  if (w1 >= (int)fAlgView.fNWires) { w1 = fAlgView.fNWires - 1; }
  
  float adc, max_adc = 0;
  for (int w = w0; w <= w1; ++w) {
    auto const* col = fAlgView.fWireDriftData[w].data();
    for (int d = d0; d <= d1; ++d) {
      adc = col[d];
      if (adc > max_adc) { max_adc = adc; }
    }
  }
  
  return max_adc;
}
// ------------------------------------------------------

//float img::DataProviderAlg::poolSum(int wire, int drift, size_t r) const
//{
//    size_t rw = r, rd = r;
//    if (!fDownscaleFullView) { rd *= fDriftWindow; }
//
//    size_t didx = getDriftIndex(drift);
//    int d0 = didx - rd; if (d0 < 0) { d0 = 0; }
//    int d1 = didx + rd; if (d1 >= (int)fNCachedDrifts) { d1 = fNCachedDrifts - 1; }
//
//    int w0 = wire - rw; if (w0 < 0) { w0 = 0; }
//    int w1 = wire + rw; if (w1 >= (int)fNWires) { w1 = fNWires - 1; }
//
//    float sum = 0;
//    for (int w = w0; w <= w1; ++w)
//    {
//        auto const * col = fWireDriftData[w].data();
//        for (int d = d0; d <= d1; ++d) { sum += col[d]; }
//    }
//
//    return sum;
//}
// ------------------------------------------------------

 std::vector<float>
img::DataProviderAlg::downscaleMax(std::size_t dst_size,
				   std::vector<float> const& adc,
				   size_t tick0) const
{
  if (debug_level>=5) std::cout << wspaces(10) << "DataProviderAlg::downscaleMax" << std::endl;   
  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;
    
    float max_adc = adc[k0] * fAlgView.fLifetimeCorrFactors[k0 + tick0];
    for (size_t k = k0 + 1; k < k1; ++k) {
      float ak = adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      if (ak > max_adc) max_adc = ak;
    }
    result[i] = max_adc;
  }
  scaleAdcSamples(result);
  return result;
}

std::vector<float>
img::DataProviderAlg::downscaleMaxMean(std::size_t dst_size,
				       std::vector<float> const& adc,
				       size_t tick0) const
{
  if (debug_level>=5) std::cout << wspaces(10) << "DataProviderAlg::downscaleMaxMean" << std::endl;   
  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;
    size_t max_idx = k0;
    float max_adc = adc[k0] * fAlgView.fLifetimeCorrFactors[k0 + tick0];
    for (size_t k = k0 + 1; k < k1; ++k) {
      float ak = adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      if (ak > max_adc) {
        max_adc = ak;
        max_idx = k;
      }
    }
    
    size_t n = 1;
    if (max_idx > 0) {
      max_adc += adc[max_idx - 1] * fAlgView.fLifetimeCorrFactors[max_idx - 1 + tick0];
      n++;
    }
    if (max_idx + 1 < adc.size()) {
      max_adc += adc[max_idx + 1] * fAlgView.fLifetimeCorrFactors[max_idx + 1 + tick0];
      n++;
    }
    
    result[i] = max_adc / n;
  }

  scaleAdcSamples(result);
  return result;
}

std::vector<float>
img::DataProviderAlg::downscaleMean(std::size_t dst_size,
				    std::vector<float> const& adc,
				    size_t tick0) const
{
  if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::downscaleMean. adc.size(), tick0 = " << adc.size() << " " << tick0 << std::endl;   
  size_t count=0; //anselmo

  size_t kStop = dst_size;
  std::vector<float> result(dst_size);
  if (adc.size() < kStop) { kStop = adc.size(); }
  for (size_t i = 0, k0 = 0; i < kStop; ++i, k0 += fDriftWindow) {
    size_t k1 = k0 + fDriftWindow;
    if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::downscaleMean, i = " << i << std::endl;   
    float sum_adc = 0;
    for (size_t k = k0; k < k1; ++k) {
      if (debug_level>=9) std::cout << wspaces(18) << "DataProviderAlg::downscaleMean, k = " << k << std::endl;
      if (k + tick0 < fAlgView.fLifetimeCorrFactors.size())
        sum_adc += adc[k] * fAlgView.fLifetimeCorrFactors[k + tick0];
      if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::downscaleMean. adc["<< k << "] = " << adc[k] << std::endl;   

    }
    result[i] = sum_adc * fDriftWindowInv;
    
    if (result[i]>0){
      if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::downscaleMean. >0 result["<< i << "] = " << result[i] << std::endl;   
      count++;
    }
  }
  if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::downscaleMean. result.size(), >0 results = " << result.size() << " " << count << std::endl;   
  scaleAdcSamples(result);

  count=0; // anselmo
  if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::downscaleMean. results after downscaling. result.size() = " << result.size() << std::endl;   
  for (size_t i=0;i<result.size();i++){
    if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::downscaleMean. result["<<i<<"] = " << result[i] << std::endl;   
    if (result[i]>0) count++;
  }
  if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::downscaleMean. results after downscaling. >0 results = " << count << std::endl;   
    
  return result;
}

//std::optional<std::vector<float>>  // anselmo
std::vector<float>
img::DataProviderAlg::setWireData(std::vector<float> const& adc, size_t wireIdx) const
{
  if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireData" << std::endl;   
  
  if (wireIdx >= fAlgView.fWireDriftData.size()) return std::vector<float>();//std::nullopt;
  auto& wData = fAlgView.fWireDriftData[wireIdx];
  
  if (fDownscaleFullView) {
    if (!adc.empty()) { return downscale(wData.size(), adc, 0); }
    else {
      return std::vector<float>();//std::nullopt;
    }
  }
  else {
    if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireData. adc.size(). fWireDriftData[wireIdx].size() = " << adc.size() << " " << wData.size()<< std::endl;             
    if (adc.empty()) { std::vector<float>();}//std::nullopt; }
    else if (adc.size() <= wData.size()){
      if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireData. 1" << std::endl;             
      return adc;
    }
    else {
      if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireData. 2" << std::endl;   
      return std::vector<float>(adc.begin(), adc.begin() + wData.size());
    }
  }
  //  return std::make_optional(wData);
  return wData;
}
// ------------------------------------------------------

bool
img::DataProviderAlg::setWireDriftData(detinfo::DetectorClocksData const& clock_data,
				       detinfo::DetectorPropertiesData const& det_prop,
				       const std::vector<my_recob::Wire>& wires,
				       unsigned int plane,
				       unsigned int tpc,
				       unsigned int cryo)
{

  if (debug_level>=1) std::cout << wspaces(2) << "DataProviderAlg::setWireDriftData. Loop over #wires = " << wires.size() << std::endl;   
  
  //  mf::LogInfo("DataProviderAlg") << "Create image for cryo:" << cryo << " tpc:" << tpc
  //				 << " plane:" << plane;
  std::cout  << "Create image for cryo:" << cryo << " tpc:" << tpc
             << " plane:" << plane << std::endl;
  
  fCryo = cryo;
  fTPC = tpc;
  fPlane = plane;
  
  fAdcSumOverThr = 0;
  fAdcAreaOverThr = 0;
  
  size_t nwires = 480;//fGeometry->Nwires(plane, tpc, cryo);   anselmo
  size_t ndrifts = 6000;//det_prop.NumberTimeSamples();    anselmo
  
  fAlgView = resizeView(clock_data, det_prop, nwires, ndrifts);

  //    auto const& channelStatus =
  //    art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();  
  const channel::ChannelStatusProvider channelStatus;
  bool allWrong = true;
  for (auto const& wire : wires) {  // loop over my_recob::Wire
    auto wireChannelNumber = wire.Channel();

    if (debug_level>=2) std::cout << wspaces(4) << "DataProviderAlg::setWireDriftData. Wire channel number = " << wireChannelNumber << std::endl;   
    if (!channelStatus.IsGood(wireChannelNumber)) { continue; }
    size_t w_idx = 0;

    if (debug_level>=2) std::cout << wspaces(4) << "DataProviderAlg::setWireDriftData. Loop over the wireIDs corresponding to that channel. wireIDs.size() = " << ChannelToWire(wireChannelNumber).size()<< std::endl;   

    //    for (auto const& id : fGeometry->ChannelToWire(wireChannelNumber)) {  // anselmo ChannelToWire returns std::vector<geo::WireID>
    for (auto const& id : ChannelToWire(wireChannelNumber)) { 
      if ((id.Plane == plane) && (id.TPC == tpc) && (id.Cryostat == cryo)) {
        w_idx = id.Wire;
        auto adc = wire.Signal();
        if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftData. Found wireID in required TPC,cryo and plane. w_idx, adc.size() = " << w_idx << " " << adc.size()<< std::endl;   
        if (adc.size() < ndrifts) {
          //mf::LogWarning("DataProviderAlg") << "Wire ADC vector size lower than NumberTimeSamples.";
          std::cout << "Wire ADC vector size lower than NumberTimeSamples." << std::endl;
          continue; // not critical, maybe other wires are OK, so continue
        }
        auto wire_data = setWireData(adc, w_idx);

        size_t l=0;
        size_t non0_values=0;;
        for (auto v : wire_data) {
          if (v!=0){
            if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::setWireDriftData. Non 0 adc values in that channel. l, v = " << l << " " << v << std::endl;   
            non0_values++;
          }
          l++;
        }

        if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftData. wire_data.size() = " << wire_data.size() << std::endl;
        if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftData. #non 0 adc values = " << non0_values << std::endl;   
        //    if (!wire_data) {
        if (wire_data.empty()) {  // anselmo
          //mf::LogWarning("DataProviderAlg") << "Wire data not set.";
          std::cout << "Wire data not set." << std::endl;
          continue; // also not critical, try to set other wires
        }
        //    fAlgView.fWireDriftData[w_idx] = *wire_data;
        fAlgView.fWireDriftData[w_idx] = wire_data;  // anselmo
        for (auto v : adc) {
          if (v >= fAdcSumThr) {
            fAdcSumOverThr += v;
            fAdcAreaOverThr++;
          }
        }
        if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftData. fAdcSumOverThr, fAdcSumAreaThr = " << fAdcSumOverThr << " " << fAdcAreaOverThr << std::endl;   
        
        fAlgView.fWireChannels[w_idx] = wireChannelNumber;
        allWrong = false;
      }
    }
  }



  if (allWrong) {
    //mf::LogError("DataProviderAlg")
    std::cout  << "Wires data not set in the cryo:" << cryo << " tpc:" << tpc << " plane:" << plane << std::endl;
    return false;
  }

  //applyBlur();
  //addWhiteNoise();
  //addCoherentNoise();
    
  return true;
}
// ------------------------------------------------------


// implemented by Anselmo
bool
img::DataProviderAlg::setWireDriftDataFromHit(const AnaHitPD& hit)
{

  if (debug_level>=1) std::cout << wspaces(2) << "DataProviderAlg::setWireDriftDataFromHit"<< std::endl;   

  static bool first=true;

  detinfo::DetectorClocksData clock_data;
  detinfo::DetectorPropertiesData det_prop;
  
  fCryo  = hit.fWireID.Cryostat;
  fTPC   = hit.fWireID.TPC;
  fPlane = hit.fWireID.Plane;

  
  fAdcSumOverThr = 0;
  fAdcAreaOverThr = 0;
  
  size_t nwires = 480;//fGeometry->Nwires(plane, tpc, cryo);   anselmo
  size_t ndrifts = 6000;//det_prop.NumberTimeSamples();    anselmo

  if (first){
    fAlgView = resizeView(clock_data, det_prop, nwires, ndrifts);
    first=false;
  }
    
  size_t w_idx = hit.fWireID.Wire;  
  //  auto adc = hit.Signal();


  size_t fPatchSizeD = 48;
  std::vector<Float_t> adc(ndrifts,0);
  for (size_t i=0;i<hit.Signal().size();i++)
    adc[i+hit.PeakTime()-fPatchSizeD/2.] = hit.Signal()[i];
  
  
  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::setWireDriftDataFromHit. Found wireID in required TPC,cryo and plane. w_idx, adc.size() = " << w_idx << " " << adc.size()<< std::endl;   
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
  
  fAlgView.fWireChannels[w_idx] = hit.fChannel;

  if (debug_level>=10) 
    for (auto v : wireData(w_idx)) {
      std::cout << wspaces(8) << "DataProviderAlg::setWireDriftDataFromHit. fAlgView.fWireDriftData[w_idx] = " << v << std::endl;   
    }

  
  //applyBlur();
  //addWhiteNoise();
  //addCoherentNoise();
    
  return true;
}
// ------------------------------------------------------


float
img::DataProviderAlg::scaleAdcSample(float val) const
{

  if (debug_level>=6) std::cout << wspaces(12) << "DataProviderAlg::scaleAdcSample" << std::endl;   

  if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::scaleAdcSample. val, Plane, fAmplCalibConst[fPlane] = " << val << " " << fPlane << " " << fAmplCalibConst[fPlane] << std::endl;   
  
  val *= fAmplCalibConst[fPlane]; // prescale by plane-to-plane calibration factors
  
  if (val < fAdcMin) { val = fAdcMin; } // saturate
  else if (val > fAdcMax) {
    val = fAdcMax;
  }
  
  return fAdcOffset +
    fAdcScale *
    (val - fAdcMin); // shift and scale to the output range, shift to the output min
}
// ------------------------------------------------------

void
img::DataProviderAlg::scaleAdcSamples(std::vector<float>& values) const
{
  if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::scaleAdcSamples. fPlane = " << fPlane << std::endl;   
  float calib = fAmplCalibConst[fPlane];
  auto* data = values.data();
  size_t k = 0, size4 = values.size() >> 2, size = values.size();
  if (debug_level>=8) std::cout << wspaces(16) << "DataProviderAlg::scaleAdcSamples. size4, size, calib = " << size4 << " " << size << " " << calib << std::endl;   
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
      if (debug_level>=9) std::cout << wspaces(18) << "DataProviderAlg::scaleAdcSamples. i, k, data[k], values[k] = " << i << " " << k << " " << data[k] << " " << values[k] << std::endl;   
    }
  while (k < size) {
    if (debug_level>=9) std::cout << wspaces(18) << "DataProviderAlg::scaleAdcSamples. k, data[k], values[k] = " << k << " " << data[k] << " " << values[k] << std::endl;   
    data[k] = scaleAdcSample(data[k]);
    ++k;
  }// do the tail
}
// ------------------------------------------------------

void
img::DataProviderAlg::applyBlur()
{
  if (fBlurKernel.size() < 2) return;
  
  size_t margin_left = (fBlurKernel.size() - 1) >> 1,
    margin_right = fBlurKernel.size() - margin_left - 1;
  
  std::vector<std::vector<float>> src(fAlgView.fWireDriftData.size());
  for (size_t w = 0; w < fAlgView.fWireDriftData.size(); ++w) {
    src[w] = fAlgView.fWireDriftData[w];
  }
    
  for (size_t w = margin_left; w < fAlgView.fWireDriftData.size() - margin_right; ++w) {
    for (size_t d = 0; d < fAlgView.fWireDriftData[w].size(); ++d) {
      float sum = 0;
      for (size_t i = 0; i < fBlurKernel.size(); ++i) {
	sum += fBlurKernel[i] * src[w + i - margin_left][d];
      }
      fAlgView.fWireDriftData[w][d] = sum;
    }
  }
}
// ------------------------------------------------------

// MUST give the same result as get_patch() in scripts/utils.py
bool
img::DataProviderAlg::patchFromDownsampledView(size_t wire,
					       float drift,
					       size_t size_w,
					       size_t size_d,
					       std::vector<std::vector<float>>& patch) const
{
  std::cout << "DataProviderAlg::patchFromDownsampledView" << std::endl;   
  int halfSizeW = size_w / 2;
  int halfSizeD = size_d / 2;
  
  int w0 = wire - halfSizeW;
  int w1 = wire + halfSizeW;
  
  size_t sd = (size_t)(drift / fDriftWindow);
  int d0 = sd - halfSizeD;
  int d1 = sd + halfSizeD;
  
  int wsize = fAlgView.fWireDriftData.size();
  for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch) {
    auto& dst = patch[wpatch];
    if ((w >= 0) && (w < wsize)) {
      auto& src = fAlgView.fWireDriftData[w];
      int dsize = src.size();
      for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch) {
	if ((d >= 0) && (d < dsize)) { dst[dpatch] = src[d]; }
	else {
	  dst[dpatch] = fAdcZero;
	}
	}
    }
    else {
      std::fill(dst.begin(), dst.end(), fAdcZero);
    }
  }
  
  return true;
}

bool
img::DataProviderAlg::patchFromOriginalView(size_t wire,
					    float drift,
					    size_t size_w,
					    size_t size_d,
					    std::vector<std::vector<float>>& patch) const
{
  if (debug_level>=3) std::cout << wspaces(6) << "DataProviderAlg::patchFromOriginalView. wire, drift = " << wire << " " << drift << std::endl;   
  int dsize = fDriftWindow * size_d;
  int halfSizeW = size_w / 2;
  int halfSizeD = dsize / 2;
  
  int w0 = wire - halfSizeW;
  int w1 = wire + halfSizeW;
  
  int d0 = int(drift) - halfSizeD;
  int d1 = int(drift) + halfSizeD;
  
  if (d0 < 0) d0 = 0;

  
  if (debug_level>=4) std::cout << wspaces(8) << "DataProviderAlg::patchFromOriginalView. w0, w1, d0, d1 = " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;   
  
  std::vector<float> tmp(dsize);
  int wsize = fAlgView.fWireDriftData.size();
  for (int w = w0, wpatch = 0; w < w1; ++w, ++wpatch) {
    if (debug_level>=5) std::cout << wspaces(10) << "DataProviderAlg::patchFromOriginalView. w, wpatch = " << w << " " << wpatch << std::endl;   
    size_t count_gt0=0;  // anselmo
    if ((w >= 0) && (w < wsize)) {
      auto& src = fAlgView.fWireDriftData[w];
      if (debug_level>=6) std::cout << wspaces(12) << "DataProviderAlg::patchFromOriginalView. fWireDriftData[w].size() = " << fAlgView.fWireDriftData[w].size() << std::endl;   
      int src_size = src.size();
      for (int d = d0, dpatch = 0; d < d1; ++d, ++dpatch) {
        if ((d >= 0) && (d < src_size)) { tmp[dpatch] = src[d]; }
        else {
          tmp[dpatch] = fAdcZero;
        }
        if (debug_level>=7) std::cout << wspaces(14) << "DataProviderAlg::patchFromOriginalView. d,dpatch,src[d]  = " << d << " " << dpatch << " " << src[d] << std::endl;   
        if (src[d]>0) count_gt0++;  // anselmo
      }
    }
    else {
      std::fill(tmp.begin(), tmp.end(), fAdcZero);
    }
    patch[wpatch] = downscale(patch[wpatch].size(), tmp, d0);
    if (debug_level>=6) std::cout << wspaces(12) << "DataProviderAlg::patchFromOriginalView. tmp.size, patch[wpatch].size(), #>0 = " << tmp.size() << " " << patch[wpatch].size() << " " << count_gt0 << std::endl;   
  }
  if (debug_level>=5) std::cout << wspaces(10) << "DataProviderAlg::patchFromOriginalView. patch.size() = " << patch.size() << std::endl;   
  return true;
}
// ------------------------------------------------------
/*
void
img::DataProviderAlg::addWhiteNoise()
{
  if (fNoiseSigma == 0) return;
  
  double effectiveSigma = scaleAdcSample(fNoiseSigma);
  if (fDownscaleFullView) effectiveSigma /= fDriftWindow;
  
  CLHEP::RandGauss gauss(fRndEngine);
  std::vector<double> noise(fAlgView.fNCachedDrifts);
  for (auto& wire : fAlgView.fWireDriftData) {
    gauss.fireArray(fAlgView.fNCachedDrifts, noise.data(), 0., effectiveSigma);
    for (size_t d = 0; d < wire.size(); ++d) {
      wire[d] += noise[d];
    }
  }
}

// ------------------------------------------------------

void
img::DataProviderAlg::addCoherentNoise()
{
  if (fCoherentSigma == 0) return;
  
  double effectiveSigma = scaleAdcSample(fCoherentSigma);
  if (fDownscaleFullView) effectiveSigma /= fDriftWindow;
  
  CLHEP::RandGauss gauss(fRndEngine);
  std::vector<double> amps1(fAlgView.fWireDriftData.size());
  std::vector<double> amps2(1 + (fAlgView.fWireDriftData.size() / 32));
  gauss.fireArray(amps1.size(), amps1.data(), 1., 0.1); // 10% wire-wire ampl. variation
  gauss.fireArray(amps2.size(), amps2.data(), 1., 0.1); // 10% group-group ampl. variation
  
  double group_amp = 1.0;
  std::vector<double> noise(fAlgView.fNCachedDrifts);
  for (size_t w = 0; w < fAlgView.fWireDriftData.size(); ++w) {
    if ((w & 31) == 0) {
      group_amp = amps2[w >> 5]; // div by 32
      gauss.fireArray(fAlgView.fNCachedDrifts, noise.data(), 0., effectiveSigma);
    } // every 32 wires
    
    auto& wire = fAlgView.fWireDriftData[w];
    for (size_t d = 0; d < wire.size(); ++d) {
      wire[d] += group_amp * amps1[w] * noise[d];
    }
  }
  }*/
// ------------------------------------------------------


// implemented by anselmo
std::vector<AnaWireID> img::DataProviderAlg::ChannelToWire(UInt_t wireChannelNumber){

  std::vector<AnaWireID> wireIDs;

  AnaWireID wireID;
  wireID.Wire=wireChannelNumber;
  wireID.Plane=2;
  wireIDs.push_back(wireID);

  return wireIDs;
}
