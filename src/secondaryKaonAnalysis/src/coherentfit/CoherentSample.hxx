#ifndef CoherentSample_hxx
#define CoherentSample_hxx

#include "pdDataClasses.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TRatioPlot.h"
#include "TRandom3.h"

class CoherentSample{
public:

  enum SampleTypeEnum {kUnassigned=0,
		       kSignalPlusBackground,
		       kSignal,
		       kSemiSignal, 
		       kBackground,
		       kSemiBackground,
		       kTrueSignal,
		       kTrueSemiSignal,
		       kTrueBackground,
		       kTrueSemiBackground};
  
  enum BackgroundModelEnum {kFitUnassigned=0,
			    kAllParFree,
			    kShift,
			    k3Par,
			    kQuadraticWidths};
    
  CoherentSample();
  CoherentSample(SampleTypeEnum type);
  virtual ~CoherentSample();

  /// Clone this object.
  CoherentSample* Clone() {
    return new CoherentSample(*this);
  }

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  CoherentSample(const CoherentSample &c);
  
public:
  
  void ComputeIntegral();
  void NormalizeHistograms();
  void NormalizeHistograms(std::vector<double> Integral);
  void UndoHistogramsNormalization() {for(int i = 0; i < (int)fh.size(); i++)fh[i]->Scale(fIntegral[i]);}
  void ChangeHistogramsToVariableBinning(const int min = 3);
  void SetCFitParameters(const CoherentSample* sample);
  void SetCFitParametersWithVariations(const CoherentSample* sample, TRandom3* r);
  std::pair<double,double> ApplyVariation(const std::pair<double,double> p, TRandom3* r);
  
  void SequentialCoherentFit(bool minos = true);
  void IncoherentFit();
  void IncoherentFitSignal();
  void IncoherentFitBackground();
  void GetInitialParValuesForCoherentFit();
  void GetInitialParValuesForSignal();
  void GetInitialParValuesForBackground();
  void CoherentFit(bool minos = true);

  void CoherentFitSignal();
  void CoherentFitOnlySignal();
  void CoherentFitBackground();
  void CoherentFitSignalPlusBackground(bool minos = true, double step = 0.001);
  void CoherentFitSignalPlusBackgroundToy(bool &result);
  void TrySoloPeaks(bool &result);
  void EstimateWhatIsLeft();
  
  void StoreCoherentFits();
  void StoreCoherentFitsSignal();
  void StoreCoherentFitsTrueSemiSignal();
  void StoreCoherentFitsBackground();
  void StoreCoherentFitsSemiBackground();

  void WriteToRootFile(const std::string& filename);
  void ReadFromRootFile(const std::string& filename);
  void ResetAllButHistograms();
  
  static void fcnSignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnOnlySignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnEstimateWhatIsLeft(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  
  CoherentSample* GetSignal() const {return fSignal;}
  void SetSignal(CoherentSample* Signal){fSignal = Signal;}
  CoherentSample* GetSemiSignal() const {return fSemiSignal;}
  void SetSemiSignal(CoherentSample* SemiSignal){fSemiSignal = SemiSignal;}
  CoherentSample* GetBackground() const {return fBackground;}
  void SetBackground(CoherentSample* Background){fBackground = Background;}
  CoherentSample* GetSemiBackground() const {return fSemiBackground;}
  void SetSemiBackground(CoherentSample* SemiBackground){fSemiBackground = SemiBackground;}
  CoherentSample* GetTrueSignal() const {return fTrueSignal;}
  void SetTrueSignal(CoherentSample* Signal){fTrueSignal = Signal;}
  CoherentSample* GetTrueSemiSignal() const {return fTrueSemiSignal;}
  void SetTrueSemiSignal(CoherentSample* SemiSignal){fTrueSemiSignal = SemiSignal;}
  CoherentSample* GetTrueBackground() const {return fTrueBackground;}
  void SetTrueBackground(CoherentSample* Background){fTrueBackground = Background;}
  CoherentSample* GetTrueSemiBackground() const {return fTrueSemiBackground;}
  void SetTrueSemiBackground(CoherentSample* SemiBackground){fTrueSemiBackground = SemiBackground;}

  SampleTypeEnum GetSampleType() const {return fType;}
  void SetSampleType(SampleTypeEnum type) {fType = type;}
  BackgroundModelEnum GetBackgroundModel() const {return fBackgroundModel;}
  void SetBackgroundModel(BackgroundModelEnum f) {fBackgroundModel = f;}

  void SetChi2Cut(const Double_t Chi2Cut) {fChi2Cut = Chi2Cut;}
  Double_t GetChi2Cut() const {return fChi2Cut;}
  
  int GetSize() const {return (int)fh.size();}
  
  std::vector<TH1F*> GetHistVector() const {return fh;}
  void AddToHistVector(TH1F* h){fh.push_back(h);}
  void ResetHistVector(){fh.clear();}

  std::vector<std::vector<TH1F*>> GetSystHistVector() const {return fSystHist;}
  
  std::vector<TF1*> GetIFitVector() const {return fIFit;}
  void AddToIFitVector(TF1* f){fIFit.push_back(f);}
  void ResetIFitVector(){fIFit.clear();}

  std::vector<TF1*> GetCFitVector() const {return fCFit;}
  void AddToCFitVector(TF1* f){fCFit.push_back(f);}
  void ResetCFitVector(){fCFit.clear();}

  TF1* GetClwFit() const {return fClwFit;}
  void SetClwFit(TF1* f){fClwFit = f;}
  TF1* GetCmpvFit() const {return fCmpvFit;}
  void SetCmpvFit(TF1* f){fCmpvFit = f;}
  TF1* GetCgwFit() const {return fCgwFit;}
  void SetCgwFit(TF1* f){fCgwFit = f;}
  
  std::vector<std::pair<double,double>> GetRRVector() const {return fRR;}
  std::vector<double> GetRRVectorValue() const;
  std::vector<double> GetRRVectorError() const;
  void AddToRRVector(std::pair<double,double> p){fRR.push_back(p);}
  void SetRRVector(std::vector<std::pair<double,double>> RR){fRR = RR;}
  void ResetRRVector(){fRR.clear();}

  std::vector<double> GetIntegralVector() const {return fIntegral;}
  void AddToIntegralVector(double Integral){fIntegral.push_back(Integral);}
  void SetIntegralVector(std::vector<double> Integral){fIntegral = Integral;}
  void ResetIntegralVector(){fIntegral.clear();}

  std::vector<double> GetIIntegralVector() const {return fIIntegral;}
  std::vector<double> GetCIntegralVector() const {return fCIntegral;}

  std::vector<std::pair<double,double>> GetIlwVector() const {return fIlw;}
  std::vector<double> GetIlwVectorValue() const;
  std::vector<double> GetIlwVectorError() const;
  void AddToIlwVector(std::pair<double,double> lw){fIlw.push_back(lw);}
  void ResetIlwVector(){fIlw.clear();}

  std::vector<std::pair<double,double>> GetImpvVector() const {return fImpv;}
  std::vector<double> GetImpvVectorValue() const;
  std::vector<double> GetImpvVectorError() const;
  void AddToImpvVector(std::pair<double,double> mpv){fImpv.push_back(mpv);}
  void ResetImpvVector(){fImpv.clear();}  

  std::vector<std::pair<double,double>> GetInormVector() const {return fInorm;}
  std::vector<double> GetInormVectorValue() const;
  std::vector<double> GetInormVectorError() const;
  void AddToInormVector(std::pair<double,double> norm){fInorm.push_back(norm);}
  void ResetInormVector(){fInorm.clear();}

  std::vector<std::pair<double,double>> GetIgwVector() const {return fIgw;}
  std::vector<double> GetIgwVectorValue() const;
  std::vector<double> GetIgwVectorError() const;
  void AddToIgwVector(std::pair<double,double> gw){fIgw.push_back(gw);}
  void ResetIgwVector(){fIgw.clear();}

  std::vector<std::pair<double,double>> GetIfwVector() const {return fIfw;}
  std::vector<double> GetIfwVectorValue() const;
  std::vector<double> GetIfwVectorError() const;
  void AddToIfwVector(std::pair<double,double> fw){fIfw.push_back(fw);}
  void ResetIfwVector(){fIfw.clear();}

  std::vector<double> GetIalphaVector() const;

  TGraphErrors* GetImpvGraph() const;
  TGraphErrors* GetIlwGraph() const;
  TGraphErrors* GetIgwGraph() const;
  TGraphErrors* GetInormGraph() const;
  TGraphErrors* GetIfwGraph() const;
  TGraphErrors* GetIconstGraph() const;
    
  std::pair<double,double> GetClwA() const {return fClwA;}
  void SetClwA(std::pair<double,double> p){fClwA = p;} 
  std::pair<double,double> GetClwB() const {return fClwB;}
  void SetClwB(std::pair<double,double> p){fClwB = p;}
  std::pair<double,double> GetClwC() const {return fClwC;}
  void SetClwC(std::pair<double,double> p){fClwC = p;}

  std::pair<double,double> GetCmpvA() const {return fCmpvA;}
  void SetCmpvA(std::pair<double,double> p){fCmpvA = p;}
  std::pair<double,double> GetCmpvB() const {return fCmpvB;}
  void SetCmpvB(std::pair<double,double> p){fCmpvB = p;}
  std::pair<double,double> GetCmpvC() const {return fCmpvC;}
  void SetCmpvC(std::pair<double,double> p){fCmpvC = p;}
  std::pair<double,double> GetCmpvD() const {return fCmpvD;}
  void SetCmpvD(std::pair<double,double> p){fCmpvD = p;}
  std::pair<double,double> GetCmpvR() const {return fCmpvR;}
  void SetCmpvR(std::pair<double,double> p){fCmpvR = p;}

  std::pair<double,double> GetCgwA() const {return fCgwA;}
  void SetCgwA(std::pair<double,double> p){fCgwA = p;}
  std::pair<double,double> GetCgwB() const {return fCgwB;}
  void SetCgwB(std::pair<double,double> p){fCgwB = p;}
  std::pair<double,double> GetCgwC() const {return fCgwC;}
  void SetCgwC(std::pair<double,double> p){fCgwC = p;}

  std::pair<double,double> GetCshift() const {return fCshift;}
  void SetCshift(std::pair<double,double> p){fCshift = p;}

  std::pair<double,double> GetClwQa() const {return fClwQa;}
  void SetClwQa(std::pair<double,double> p){fClwQa = p;}

  std::pair<double,double> GetCgwQa() const {return fCgwQa;}
  void SetCgwQa(std::pair<double,double> p){fCgwQa = p;} 
  
  std::pair<double,double> GetCnorm() const {return fCnorm;}
  void SetCnorm(std::pair<double,double> p){fCnorm = p;}

  std::pair<double,double> GetCconst() const {return fCconst;}
  void SetCconst(std::pair<double,double> p){fCconst = p;}

  void SetHistogramMarker(int n){for(int i = 0; i < (int)fh.size(); i++)fh[i]->SetMarkerStyle(n);}
  void SetHistogramColor(int n){for(int i = 0; i < (int)fh.size(); i++){fh[i]->SetMarkerColor(n);fh[i]->SetLineColor(n);}}
  void SetHistogramLineWidth(int n){for(int i = 0; i < (int)fh.size(); i++){fh[i]->SetLineWidth(n);}}
  void UseCurrentStyle(){for(int i = 0; i < (int)fh.size(); i++){fh[i]->UseCurrentStyle();}}
  
  void SetIFitStyle(int n){for(int i = 0; i < (int)fIFit.size(); i++)fIFit[i]->SetLineStyle(n);}
  void SetIFitColor(int n){for(int i = 0; i < (int)fIFit.size(); i++)fIFit[i]->SetLineColor(n);}

  void SetCFitStyle(int n){for(int i = 0; i < (int)fCFit.size(); i++)fCFit[i]->SetLineStyle(n);}
  void SetCFitColor(int n){for(int i = 0; i < (int)fCFit.size(); i++)fCFit[i]->SetLineColor(n);}

  TGraphErrors* GetMPVErrorBandOld();
  TGraphAsymmErrors* GetMPVErrorBand(const double min = 1, const double max = 31, const double step = 0.01);

  TRatioPlot* GetFitRatio(const int ibin);

  void CopyHistogramsBinning(std::vector<TH1F*> vh);

  std::vector<double> GetLikelihoodVector() const {return fPartialLikelihood;}

  double GetCorrelationMatrixElement(int i, int j){if(i < 0 || i > 17 || j < 0 || j > 17)return -999.;
                                                   else return fcorr_matrix[i][j];}
  
private:

  CoherentSample* fSignal;
  CoherentSample* fSemiSignal;
  CoherentSample* fBackground;
  CoherentSample* fSemiBackground;
  CoherentSample* fTrueSignal;
  CoherentSample* fTrueSemiSignal;
  CoherentSample* fTrueBackground;
  CoherentSample* fTrueSemiBackground;
  
  SampleTypeEnum      fType;
  BackgroundModelEnum fBackgroundModel;

  Double_t fChi2Cut;
  
  TMinuit* fMinuit;
  
  std::vector<TH1F*> fh;
  std::vector<std::vector<TH1F*>> fSystHist;

  std::vector<double> fIntegral;
  std::vector<double> fIIntegral;
  std::vector<double> fCIntegral;

  std::vector<TF1*> fIFit;
  std::vector<TF1*> fCFit;

  TF1* fClwFit;
  TF1* fCmpvFit;
  TF1* fCgwFit;

  std::vector<std::pair<double,double>> fRR;

  //incoherent results
  std::vector<std::pair<double,double>> fIlw;
  std::vector<std::pair<double,double>> fImpv;
  std::vector<std::pair<double,double>> fInorm;
  std::vector<std::pair<double,double>> fIgw;
  std::vector<std::pair<double,double>> fIfw;
  std::vector<std::pair<double,double>> fIalpha;
  std::vector<std::pair<double,double>> fIconst;

  //coherent signal results
  std::pair<double,double> fClwA;
  std::pair<double,double> fClwB;
  std::pair<double,double> fClwC;
  
  std::pair<double,double> fCmpvA;
  std::pair<double,double> fCmpvB;
  std::pair<double,double> fCmpvC;
  std::pair<double,double> fCmpvD;
  std::pair<double,double> fCmpvR;

  std::pair<double,double> fCgwA;
  std::pair<double,double> fCgwB;
  std::pair<double,double> fCgwC;

  std::pair<double,double> fCshift;
  std::pair<double,double> fClwQa;
  std::pair<double,double> fCgwQa;

  std::pair<double,double> fCnorm;

  std::pair<double,double> fCconst;

  double fCmpvA_error_neg;
  double fCmpvA_error_pos;
  double fCmpvB_error_neg;
  double fCmpvB_error_pos;
  double fCmpvC_error_neg;
  double fCmpvC_error_pos;

  double fLikelihood;
  std::vector<double> fPartialLikelihood;

  double fcorr_matrix[18][18];
};

#endif
