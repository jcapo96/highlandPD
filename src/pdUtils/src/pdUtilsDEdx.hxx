#ifndef pdUtilsDEdx_h
#define pdUtilsDEdx_h

#include "pdDataClasses.hxx"
#include <utility>
#include <vector>

class TGraph;
class TMultiGraph;
class TH1F;

namespace pdAnaUtils {

  struct DEdxFreeRangeFitResult {
    Float_t logLikelihood = -999.f;
    Float_t bestOffsetCm = -999.f;
    Float_t momentumGeV = -999.f;
    Float_t meanDedxBias = -999.f;   // Gaussian mean of (measured - Bethe-Bloch) dEdx distribution [MeV/cm]
    Float_t sigmaDedxBias = -999.f;  // Gaussian sigma of (measured - Bethe-Bloch) dEdx distribution [MeV/cm]
    Int_t dedxFitOk = -1;
  };

  Float_t GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG, int skipHitsFirst = 1, int skipHitsLast = 1,
                            double dedxMinMeVcm = 0., double dedxMaxMeVcm = 0.);
  Float_t GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR, int skipHitsFirst = 1,
                                   int skipHitsLast = 1, double dedxMinMeVcm = 0., double dedxMaxMeVcm = 0.);
  Float_t dEdxLikelihood(TGraph* tg, TGraph* tg_ke, Float_t mass);

  std::pair<Float_t, Float_t> GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG, int skipHitsFirst = 1,
                                                         int skipHitsLast = 1, double dedxMinMeVcm = 0.,
                                                         double dedxMaxMeVcm = 0., double pdfPathCm = 0.65);
  std::pair<Float_t, Float_t> GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                                int skipHitsFirst = 1, int skipHitsLast = 1,
                                                                double dedxMinMeVcm = 0., double dedxMaxMeVcm = 0.,
                                                                double pdfPathCm = 0.65);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax = 450.,
                                                       double step = 1.0, int minInteriorPoints = 15,
                                                       int skipHitsFirst = 3, int skipHitsLast = 3,
                                                       double dedxMinMeVcm = 0.5, double dedxMaxMeVcm = 5.0,
                                                       double pdfPathCm = 0.65);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                               double Lmax = 500., double step = 0.5,
                                                               int minInteriorPoints = 2, int skipHitsFirst = 1,
                                                               int skipHitsLast = 1, double dedxMinMeVcm = 0.,
                                                               double dedxMaxMeVcm = 0., double pdfPathCm = 0.65);

  TMultiGraph* MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax, double step, int skipHitsFirst,
                                                   int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm,
                                                   int minInteriorPoints, double pdfPathCm,
                                                   const char* xAxisTitle = nullptr);

  TH1F* MakePionFreeRangeDedxBiasHistogram(AnaParticlePD* part, double Lmax, double step, int skipHitsFirst,
                                           int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm,
                                           int minInteriorPoints, double pdfPathCm, const char* histTitle = nullptr);

  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke, Float_t mass);
  DEdxFreeRangeFitResult dEdxLikelihoodFreeRangeFit(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0,
                                                    double Lmax, double step, double measuredTrackLengthCm,
                                                    bool computeMomentum, double pdfPathCm = 0.65);

  /// For π (PDG 211): same free-range scan as GetdEdxLikelihoodFreeRangeFit, map each trial offset L to
  /// momentum via R_eff = measured_length + L and RangeCmToMomentumGeV. Duplicate p bins keep max log L.
  /// Returns false if fewer than two usable (p, logL) points.
  bool BuildPionFreeRangeLogLikelihoodVsMomentumCurve(AnaParticlePD* part, double Lmax, double step,
                                                       int minInteriorPoints, int skipHitsFirst, int skipHitsLast,
                                                       double dedxMinMeVcm, double dedxMaxMeVcm, double pdfPathCm,
                                                       std::vector<double>& pGeV, std::vector<double>& logL);

  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
  double KineticEnergyMeVFromResidualRangeCm(TGraph* tg_ke, double range_cm);
  double GetWmax(double KE, double mass);
  double GetLandauXi(double KE, double dx, double mass);
  double dEdxPDF(double *x, double *par);

}

#endif
