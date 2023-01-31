#ifndef CoherentFitUtils_h
#define CoherentFitUtils_h

#include "CoherentFit.hxx"
#include "CoherentFitAlgorithm.hxx"

#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"

namespace CoherentFitUtils{
  
  TH1F* CreateHistogram(double bin_min, double bin_max, double bin_width);
  
  
  TH1F* GetHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill,
				      const double RMIN, const double RMAX,
				      const double Chi2Cut,
				      double bin_min = 0, double bin_max = 50, double bin_width = 0.1);

  TH1F* GetSignalHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
					    const double RMIN, const double RMAX,
					    const double Chi2Cut);
  
  TH1F* GetBackgroundHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
						const double RMIN, const double RMAX,
						const double Chi2Cut);

  TH1F* GetSemiSignalHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
						const double RMIN, const double RMAX,
						const double Chi2Cut);
  
  TH1F* GetAllKaonsHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
					      const double RMIN, const double RMAX,
					      const double Chi2Cut);
  
  TH1F* GenerateBackgroundHistogramFromTrueSignal(TTree* t, AnaSpillB* spill, TH1F* ha,
						  const double RMIN, const double RMAX,
						  const double Chi2Cut,
						  const std::vector<double> shift);
  
  TH1F* GetHistogramFromResRangeSliceFromFlatTree(TTree* t,
						  const double RMIN, const double RMAX,
						  double bin_min = 0, double bin_max = 50, double bin_width = 0.1,
						  std::string cut = "",
						  const bool apply_toy_weights = false, const bool apply_toy_variations = false, int itoy = -1);
  
  TH1F* GetToyHistogramFromResRangeSlice(TTree* t,
					 const double RMIN, const double RMAX,
					 double bin_min = 0, double bin_max = 50, double bin_width = 0.1,
					 const bool apply_toy_weights = false, const bool apply_toy_variations = false,
					 const int itoy = -1);

  TH1F* ChangeHistogramToVariableBinning(TH1F* h_original, const double I, const int min);

  void CopyHistogramBinning(TH1F* h_original, TH1F* h_copy);

  Double_t GetFunctionNormalizationInsideHistogramBoundaries(const TH1F* h, const TF1* f);
  
  TF1* GausnFit(TH1F* h, bool use_poisson = true);
  TF1* LangausFit(TH1F* h, CoherentSample::SampleTypeEnum sample, bool use_poisson = true);

 
  Double_t Langaus(Double_t *x, Double_t *par);
  Double_t LangausPlusConstant(Double_t *x, Double_t *par);
  Double_t DoubleLangaus(Double_t *x, Double_t *par);

  Double_t ABCParametrization(Double_t *x, Double_t *par);
  Double_t QuadraticABCParametrization(Double_t *x, Double_t *par);
  Double_t ABCDRParametrization(Double_t *x, Double_t *par);
  Double_t ABCDRDerivativeA(Double_t *x, Double_t *par);
  Double_t ABCDRDerivativeB(Double_t *x, Double_t *par);
  Double_t ABCDRDerivativeC(Double_t *x, Double_t *par);
  Double_t ABCDRDerivativeD(Double_t *x, Double_t *par);
  Double_t ABCDRDerivativeR(Double_t *x, Double_t *par);

  void GetABCParametrization(double &A, double &B, double &C,
			     std::vector<std::pair<double,double>> X, std::vector<std::pair<double,double>> Y,
			     bool equal_weights = false, bool draw_plot = false);

  void GetABCDRParametrization(double &A, double &B, double &C, double &D, double &R,
			       std::vector<std::pair<double,double>> X, std::vector<std::pair<double,double>> Y,
			       bool equal_weights = false, bool draw_plot = false);

  double ComputeLikelihood(TH1F* h, TF1* f, double integral);
  double NormRegularization(const std::vector<double> norm, const std::vector<double> integral);
  
  void SeparatePair(std::vector<std::pair<double,double>> p, std::vector<double> &v1, std::vector<double> &v2);

  void DeleteOutliers(std::vector<double> &x, std::vector<double> &x_error, std::vector<double> &y, std::vector<double> &y_error);

  bool IsTrueSignal(AnaParticlePD* part);

  TGraphErrors* GetGraph(std::vector<std::pair<double,double>> x, std::vector<std::pair<double,double>> y);
}

#endif
