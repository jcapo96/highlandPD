// Class:       PointIdAlg
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),         from DUNE, CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),               from DUNE, FNAL/NCBJ, since May 2016
//              P.Plonski,                                    from DUNE, WUT,       since May 2016
//
//
// Algorithm for making 2D image-like data from recob::Wire's. Used by CNN codes for training data
// preparation and application of trained models to new data. Also used by PMA to keep images of
// 2D projections used for the track validation.
//

#ifndef DataProviderAlg_h
#define DataProviderAlg_h

// Framework includes
//#include "fhiclcpp/types/Atom.h"
//#include "fhiclcpp/types/Sequence.h"
//#include "fhiclcpp/types/Table.h"

// LArSoft includes
//#include "larcorealg/Geometry/GeometryCore.h"
//#include "lardataobj/RecoBase/Wire.h"
//#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "CalorimetryAlg.hxx"

//#include "CLHEP/Random/JamesRandom.h" // for testing on noise, not used by any reco
 
// ROOT & C++
#include <memory>
#include <set>


const int debug_level=0;

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace my_recob{
  class Wire{
  public:
    
    Wire(){}

    std::vector<Float_t> fSignal;
    UInt_t fChannel;

    UInt_t Channel() const {return fChannel;}
    std::vector<Float_t> Signal() const {return fSignal;}
    

  };


  class Track{
  public:
    
    Track(AnaParticleB& part){
      fID = part.UniqueID;
      fEnd = TVector3(part.PositionEnd[0],   part.PositionEnd[1],   part.PositionEnd[2]);
      fLoc = TVector3(part.PositionStart[0], part.PositionStart[1], part.PositionStart[2]);
      fPart = &part;
    }
    
    Int_t fID;
    TVector3 fEnd;
    TVector3 fLoc;
    AnaParticleB* fPart;
    Float_t fLength;
    
    Int_t ID(){return fID;}
    Float_t Length(){return fLength;}
    
    
  };

}

namespace channel {
  class ChannelStatusProvider{
  public:
    //{kDISCONNECTED=0, kDEAD=1, kLOWNOISE=2, kNOISY=3, kGOOD=4, kUNKNOWN=5};
    bool IsPresent(unsigned int channel) const { return fStatus == 0 ? false : true; }
    
    bool IsBad(unsigned int channel) const 
    { return fBadChannels.count(channel) > 0; }
    
    bool IsNoisy(unsigned int channel) const 
    { return fNoisyChannels.count(channel) > 0; }
    
    bool IsGood(unsigned int channel) const {
      //      return IsPresent(channel) && !IsBad(channel) && !IsNoisy(channel);  // anselmo
      return true;
    }

    ChannelStatusProvider(){};
    ~ChannelStatusProvider(){};
  
  private:
    std::set<unsigned int> fBadChannels;
    std::set<unsigned int> fNoisyChannels;
    unsigned int fStatus;
  };
  
}



namespace img {
  class DataProviderAlg;
  struct DataProviderAlgView {
    unsigned int fNWires;
    unsigned int fNDrifts;
    unsigned int fNScaledDrifts;
    unsigned int fNCachedDrifts;
    std::vector<unsigned int> fWireChannels;
    std::vector<std::vector<float>> fWireDriftData;
    std::vector<float> fLifetimeCorrFactors;
  };
}
    
class img::DataProviderAlg {
public:
  enum EDownscaleMode { kMax = 1, kMaxMean = 2, kMean = 3 };
  
  /*struct Config {
    fhicl::Table<calo::CalorimetryAlg::Config> CalorimetryAlg{
      Name("CalorimetryAlg"),
	Comment("Used to eliminate amplitude variation due to electron lifetime.")};
    
    fhicl::Atom<float> AdcMax{Name("AdcMax"), Comment("Saturation max value")};
    fhicl::Atom<float> AdcMin{Name("AdcMin"), Comment("Saturation min value")};
    fhicl::Atom<float> OutMax{Name("OutMax"), Comment("Output max value")};
    fhicl::Atom<float> OutMin{Name("OutMin"), Comment("Output min value")};
   
    fhicl::Atom<bool> CalibrateAmpl{Name("CalibrateAmpl"),
	Comment("Calibrate ADC values with CalAmpConstants")};
   
    fhicl::Atom<bool> CalibrateLifetime{Name("CalibrateLifetime"),
	Comment("Calibrate ADC values with the electron lifetime")};
    
    fhicl::Atom<unsigned int> DriftWindow{Name("DriftWindow"),
	Comment("Downsampling window (in drift ticks).")};
    
    fhicl::Atom<std::string> DownscaleFn{Name("DownscaleFn"), Comment("Downsampling function")};
    
    fhicl::Atom<bool> DownscaleFullView{
      Name("DownscaleFullView"),
	Comment("Downsample full view (faster / lower location precision)")};
  
    fhicl::Sequence<float> BlurKernel{Name("BlurKernel"), Comment("Blur kernel in wire direction")};
   
    fhicl::Atom<float> NoiseSigma{Name("NoiseSigma"), Comment("White noise sigma")};
   
    fhicl::Atom<float> CoherentSigma{Name("CoherentSigma"), Comment("Coherent noise sigma")};
    };*/
   
  /*DataProviderAlg(const fhicl::ParameterSet& pset)
    : DataProviderAlg(fhicl::Table<Config>(pset, {})())
  {}
   
  DataProviderAlg(const Config& config);*/

  DataProviderAlg();
   
  virtual ~DataProviderAlg();

  // anselmo
  std::vector<AnaWireID> ChannelToWire(UInt_t wireChannelNumber);
  
  bool setWireDriftData(const detinfo::DetectorClocksData& clock_data,
			const detinfo::DetectorPropertiesData& det_prop,
			const std::vector<my_recob::Wire>&
			wires, // once per plane: setup ADC buffer, collect & downscale ADC's
			unsigned int plane,
			unsigned int tpc,
			unsigned int cryo);

  bool setWireDriftDataFromHit(const AnaHitPD& hit);  // anselmo
  
  std::vector<float> const&
  wireData(size_t widx) const
  {
    return fAlgView.fWireDriftData[widx];
  }

  //anselmo
  std::vector<std::vector<Float_t> >& wireData()
  {
    return fAlgView.fWireDriftData;
  }


  
  std::vector<std::vector<float>>
  getPatch(size_t wire, float drift, size_t patchSizeW, size_t patchSizeD) const
  {
    bool ok = false;
    std::vector<std::vector<float>> patch;
    if (fDownscaleFullView) {
      ok = patchFromDownsampledView(wire, drift, patchSizeW, patchSizeD, patch);
    }
    else {
      ok = patchFromOriginalView(wire, drift, patchSizeW, patchSizeD, patch);
    }
    
    if (ok)
      return patch;
    //    throw cet::exception("img::DataProviderAlg") << "Patch filling failed." << std::endl;
  }
  
  float
  getPixelOrZero(int wire, int drift) const
  {
    size_t didx = getDriftIndex(drift), widx = (size_t)wire; 
    
    if ((widx < fAlgView.fWireDriftData.size()) && (didx < fAlgView.fNCachedDrifts)) {
      return fAlgView.fWireDriftData[widx][didx];
    }
    return 0;
  }
  
  double
  getAdcSum() const
  {
    return fAdcSumOverThr;
  }
  size_t
  getAdcArea() const
  {
    return fAdcAreaOverThr;
  }

  float poolMax(int wire, int drift, size_t r = 0) const;

  //  float poolSum(int wire, int drift, size_t r = 0) const;

  unsigned int
  Cryo() const
  {
    return fCryo;
  }
  unsigned int
  TPC() const
  {
    return fTPC;
  }
  unsigned int
  Plane() const
  {
    return fPlane;
  }
  
  unsigned int
  NWires() const
  {
    return fAlgView.fNWires;
  }
  unsigned int
  NScaledDrifts() const
  {
    return fAlgView.fNScaledDrifts;
  }
  unsigned int
  NCachedDrifts() const
  {
    return fAlgView.fNCachedDrifts;
  }
  unsigned int
  DriftWindow() const
  {
    return fDriftWindow;
  }
  
  float
  ZeroLevel() const
  {
    return fAdcZero;
  }
  
  double
  LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
		     detinfo::DetectorPropertiesData const& det_prop,
		     double tick) const
  {
    //    return fCalorimetryAlg.LifetimeCorrection(clock_data, det_prop, tick);
    return fCalorimetryAlg.LifetimeCorrection(tick);
  }

protected:
  DataProviderAlgView fAlgView;
  EDownscaleMode fDownscaleMode;
  //std::function<void (std::vector<float> &, std::vector<float> const &, size_t)> fnDownscale;
  
  size_t fDriftWindow;
  bool fDownscaleFullView;
  float fDriftWindowInv;
  
  std::vector<float> downscaleMax(std::size_t dst_size,
				  std::vector<float> const& adc,
				  size_t tick0) const;
  std::vector<float> downscaleMaxMean(std::size_t dst_size,
				      std::vector<float> const& adc,
				      size_t tick0) const;
  std::vector<float> downscaleMean(std::size_t dst_size,
				   std::vector<float> const& adc,
				   size_t tick0) const;
  std::vector<float>
  downscale(std::size_t dst_size, std::vector<float> const& adc, size_t tick0) const
  {    
    if (debug_level>=6) std::cout << "            DataProviderAlg::downscale" << std::endl;   
    switch (fDownscaleMode) {
    case img::DataProviderAlg::kMean: return downscaleMean(dst_size, adc, tick0);
    case img::DataProviderAlg::kMaxMean: return downscaleMaxMean(dst_size, adc, tick0);
    case img::DataProviderAlg::kMax: return downscaleMax(dst_size, adc, tick0);
    }
    //throw cet::exception("img::DataProviderAlg") << "Downscale mode not supported." << std::endl;
  }

  size_t
  getDriftIndex(float drift) const
  {
    if (fDownscaleFullView)
      return (size_t)(drift * fDriftWindowInv);
    else
      return (size_t)drift;
  }
  
  //  std::optional<std::vector<float>> setWireData(std::vector<float> const& adc,   // anselmo
  std::vector<float> setWireData(std::vector<float> const& adc,
						size_t wireIdx) const;

  bool patchFromDownsampledView(size_t wire,
				float drift,
				size_t size_w,
				size_t size_d,
				std::vector<std::vector<float>>& patch) const;
  bool patchFromOriginalView(size_t wire,
			     float drift,
			     size_t size_w,
			     size_t size_d,
			     std::vector<std::vector<float>>& patch) const;
  
  virtual DataProviderAlgView resizeView(detinfo::DetectorClocksData const& clock_data,
					 detinfo::DetectorPropertiesData const& det_prop,
					 size_t wires,
					 size_t drifts);

  // Calorimetry needed to equalize ADC amplitude along drift:
  calo::CalorimetryAlg fCalorimetryAlg;
  /*   
  // Geometry and detector properties:
  geo::GeometryCore const* fGeometry;
  */
private:
  float scaleAdcSample(float val) const;
  void scaleAdcSamples(std::vector<float>& values) const;
  std::vector<float> fAmplCalibConst;
  bool fCalibrateAmpl, fCalibrateLifetime;
  unsigned int fCryo = 9999, fTPC = 9999, fPlane = 9999;
  float fAdcMax, fAdcMin, fAdcScale, fAdcOffset, fAdcZero;
  double fAdcSumOverThr, fAdcSumThr;
  size_t fAdcAreaOverThr;

  //  CLHEP::HepJamesRandom fRndEngine;

  void applyBlur();
  std::vector<float> fBlurKernel; // blur not applied if empty
    /*  
  void addWhiteNoise();
  float fNoiseSigma; // noise not added if sigma=0

  void addCoherentNoise();
  float fCoherentSigma; // noise not added if sigma=0
*/
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
 
#endif
