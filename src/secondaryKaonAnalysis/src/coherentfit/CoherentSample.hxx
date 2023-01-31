#ifndef CoherentSample_hxx
#define CoherentSample_hxx

#include "pdDataClasses.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

class CoherentSample{
public:

  enum SampleTypeEnum {kUnassigned=0,
		       kSignalPlusBackground,
		       kSignal,
		       kSemiSignal, 
		       kBackground,
		       kTrueSignal,
		       kTrueSemiSignal,
		       kTrueBackground};
  
  enum BackgroundModelEnum {kFitUnassigned=0,
			    kAllParFree,
			    kShift,
			    k3Par,
			    kQuadraticWidths};
    
  CoherentSample();
  CoherentSample(SampleTypeEnum type);
  virtual ~CoherentSample(){}

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
  void SetCFitParametersWithVariations(const CoherentSample* sample, TRandom3* r, const double sigma,
				       const bool apply_all_var = true,
				       const bool apply_only_lw_var = false, 
				       const bool apply_only_mpv_var = false, 
				       const bool apply_only_gw_var = false,
				       const bool apply_only_shift_var = false,
				       const bool apply_only_norm_var = false);
  
  void SequentialCoherentFit();
  void IncoherentFit();
  void GetInitialParValuesForCoherentFit(bool equal_weights = false, bool draw_fits = false);
  void CoherentFit();

  void CoherentFitSignal();
  void CoherentFitSignalCheb();
  void CoherentFitSignalLag();
  void CoherentFitSignalLagMPV();
  //void CoherentFitSignalAlpha();
  //void CoherentFitBackgroundAllFree();
  //void CoherentFitSignalPlusBackgroundAllFree();
  //void CoherentFitBackgroundShift();
  //void CoherentFitSignalPlusBackgroundShift();
  //void CoherentFitBackground3Par();
  void CoherentFitBackgroundQuadraticWidths();
  void CoherentFitBackgroundCheb();
  void CoherentFitBackgroundLag();
  void CoherentFitBackgroundLagMPV();
  void CoherentFitBackgroundMiau();
  void CoherentFitSignalPlusBackgroundQuadraticWidths();
  void CoherentFitSignalPlusBackgroundCheb();
  void CoherentFitSignalPlusBackgroundLag();
  void CoherentFitSignalPlusBackgroundLagMPV();
  void CoherentFitSignalPlusBackgroundMiau();
  
  void StoreCoherentFits();

  void WriteToRootFile(const std::string& filename);
  void ReadFromRootFile(const std::string& filename);
  void ResetAllButHistograms();
  
  static void fcnSignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnSignalAlpha(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnBackgroundAllFree(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnSignalPlusBackgroundAllFree(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnBackgroundShift(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnSignalPlusBackgroundShift(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //static void fcnBackground3Par(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackgroundQuadraticWidths(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackgroundQuadraticWidths(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackgroundCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackgroundLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackgroundLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnBackgroundMiau(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackgroundCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackgroundLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackgroundLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static void fcnSignalPlusBackgroundMiau(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  
  CoherentSample* GetSignal() const {return fSignal;}
  void SetSignal(CoherentSample* Signal){fSignal = Signal;}
  CoherentSample* GetBackground() const {return fBackground;}
  void SetBackground(CoherentSample* Background){fBackground = Background;}
  CoherentSample* GetTrueSignal() const {return fTrueSignal;}
  void SetTrueSignal(CoherentSample* Signal){fTrueSignal = Signal;}
  CoherentSample* GetTrueBackground() const {return fTrueBackground;}
  void SetTrueBackground(CoherentSample* Background){fTrueBackground = Background;}

  SampleTypeEnum GetSampleType() const {return fType;}
  void SetSampleType(SampleTypeEnum type) {fType = type;}
  BackgroundModelEnum GetBackgroundModel() const {return fBackgroundModel;}
  void SetBackgroundModel(BackgroundModelEnum f) {fBackgroundModel = f;}

  void SetChi2Cut(const Double_t Chi2Cut) {fChi2Cut = Chi2Cut;}
  Double_t GetChi2Cut() const {return fChi2Cut;}
  
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

  std::vector<double> GetIalphaVector() const;

  TGraphErrors* GetImpvGraph() const;
  TGraphErrors* GetIlwGraph() const;
  TGraphErrors* GetIgwGraph() const;
  TGraphErrors* GetInormGraph() const;
  TGraphErrors* GetIfwGraph() const;
    
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

  void SetHistogramMarker(int n){for(int i = 0; i < (int)fh.size(); i++)fh[i]->SetMarkerStyle(n);}
  void SetHistogramColor(int n){for(int i = 0; i < (int)fh.size(); i++){fh[i]->SetMarkerColor(n);fh[i]->SetLineColor(n);}}
  void SetHistogramLineWidth(int n){for(int i = 0; i < (int)fh.size(); i++){fh[i]->SetLineWidth(n);}}
  void UseCurrentStyle(){for(int i = 0; i < (int)fh.size(); i++){fh[i]->UseCurrentStyle();}}
  
  void SetIFitStyle(int n){for(int i = 0; i < (int)fIFit.size(); i++)fIFit[i]->SetLineStyle(n);}
  void SetIFitColor(int n){for(int i = 0; i < (int)fIFit.size(); i++)fIFit[i]->SetLineColor(n);}

  void SetCFitStyle(int n){for(int i = 0; i < (int)fCFit.size(); i++)fCFit[i]->SetLineStyle(n);}
  void SetCFitColor(int n){for(int i = 0; i < (int)fCFit.size(); i++)fCFit[i]->SetLineColor(n);}

  TGraphErrors* GetMPVErrorBand();

  void CopyHistogramsBinning(std::vector<TH1F*> vh);
  
private:

  CoherentSample* fSignal;
  CoherentSample* fBackground;
  CoherentSample* fTrueSignal;
  CoherentSample* fTrueBackground;
  
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
};

#endif
