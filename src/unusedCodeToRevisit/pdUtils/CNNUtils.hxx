#ifndef CNNUtils_h
#define CNNUtils_h


#include "pdDataClasses.hxx"
#include "pdCalorimetryUtils.hxx"
#include "tf_graph.hxx"
#include "timeUtils.hxx"

const int debug_levelA=2;

struct DataProviderAlgView {
  unsigned int fNWires;
  unsigned int fNDrifts;
  unsigned int fNScaledDrifts;
  unsigned int fNCachedDrifts;
  std::vector<unsigned int> fWireChannels;
  std::vector<std::vector<float>> fWireDriftData;
  std::vector<float> fLifetimeCorrFactors;
};

namespace nnet2{
  class ModelInterface {
  public:
    virtual ~ModelInterface() {}
    
    virtual std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) = 0;
    virtual std::vector<std::vector<float>> Run(
                                                std::vector<std::vector<std::vector<float>>> const& inps,
                                                int samples = -1);


    virtual void SetTimeUtils(timeUtils* time)=0;
  protected:
    std::string findFile(const char* fileName) const;
  };


  class TfModelInterface : public ModelInterface {
  public:
    TfModelInterface(const char* modelFileName);
    
    std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                        //                                      int samples = -1) override;
                                        int samples = -1);
    //  std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) override;
    std::vector<float> Run(std::vector<std::vector<float>> const& inp2d);

    void SetTimeUtils(timeUtils* time){g->_tutils=time;}
    
  private:
    std::unique_ptr<tf::Graph> g; // network graph
  };
  
}


class CNNUtils{
public:
  
  
  CNNUtils();

  virtual ~CNNUtils();

  enum EDownscaleMode { kMax = 1, kMaxMean = 2, kMean = 3 };
  
  //  void produce(std::vector<AnaHitPD*>& hits); //override;
  void produce(std::vector<AnaHitPD*>& hits, std::vector<AnaWireCNN>& wires);
  void FillHits(const AnaEventPD& evt, std::vector<AnaHitPD*>& hitPtrList);


  std::vector<std::vector<float>> predictIdVectors(std::vector<std::pair<unsigned int, float>> points);

  
  bool setWireDriftDataFromHit(const AnaHitPD& hit);


  bool setWireDriftData(const std::vector<AnaWireCNN>& wires, Int_t TPC, Int_t plane);
  
  
  void resizeView(size_t wires,size_t drifts);

  std::vector<float> setWireData(std::vector<float> const& adc, size_t wireIdx) const;


  std::vector<float> const& wireData(size_t widx) const{ return fAlgView.fWireDriftData[widx];}


  std::vector<float> downscale       (std::size_t dst_size, std::vector<float> const& adc,size_t tick0, size_t adc_gt0=0) const;    
  std::vector<float> downscaleMean   (std::size_t dst_size, std::vector<float> const& adc,size_t tick0, size_t adc_gt0) const;
  std::vector<float> downscaleMax    (std::size_t dst_size, std::vector<float> const& adc,size_t tick0) const{return std::vector<float>();}
  std::vector<float> downscaleMaxMean(std::size_t dst_size, std::vector<float> const& adc,size_t tick0) const{return std::vector<float>();}

                                  


  
  float scaleAdcSample(float val) const;
  void scaleAdcSamples(std::vector<float>& values, size_t adc_gt0) const;

  bool bufferPatch(size_t wire, float drift, std::vector<std::vector<float>>& patch) const;
  bool patchFromOriginalView(size_t wire,float drift,size_t size_w,size_t size_d,std::vector<std::vector<float>>& patch) const;
  bool patchFromDownsampledView(size_t wire,float drift,size_t size_w,size_t size_d,std::vector<std::vector<float>>& patch) const{return true;}

  void printTime(const std::string& m);

  void ComputeParticleCNN(AnaParticlePD& part);


  size_t GetWireIndex(Int_t wire);
  size_t GetWireTPC(Int_t wire);
  size_t GetWirePlane(Int_t wire);
  
public:
  timeUtils* _tutils;
  
private:    



  size_t fBatchSize;

  DataProviderAlgView fAlgView;
  EDownscaleMode fDownscaleMode;
  //std::function<void (std::vector<float> &, std::vector<float> const &, size_t)> fnDownscale;


  float fDriftWindowInv;
  
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


  
  pdCalorimetryUtils fCalo;



  std::string fNNetModelFilePath;
  std::vector<std::string> fNNetOutputs;
  nnet2::ModelInterface* fNNet;

  
  int fPlane=2;
};
 
  
#endif
