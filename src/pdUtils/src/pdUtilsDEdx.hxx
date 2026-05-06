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
                                                       double pdfPathCm = 0.65,
                                                       bool computeDedxBiasDiagnostics = true,
                                                       double scanStepFineCm = 0.,
                                                       double lowPMomentumRefineGeV = 0.2);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                               double Lmax = 500., double step = 0.5,
                                                               int minInteriorPoints = 2, int skipHitsFirst = 1,
                                                               int skipHitsLast = 1, double dedxMinMeVcm = 0.,
                                                               double dedxMaxMeVcm = 0., double pdfPathCm = 0.65,
                                                               bool computeDedxBiasDiagnostics = true,
                                                               double scanStepFineCm = 0.,
                                                               double lowPMomentumRefineGeV = 0.2);

  TMultiGraph* MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax, double step, int skipHitsFirst,
                                                   int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm,
                                                   int minInteriorPoints, double pdfPathCm,
                                                   const char* xAxisTitle = nullptr);
  TGraph* MakePionBetheBlochDedxVsRRReference(double rrMinCm, double rrMaxCm, double rrStepCm = 0.5);

  TH1F* MakePionFreeRangeDedxBiasHistogram(AnaParticlePD* part, double Lmax, double step, int skipHitsFirst,
                                           int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm,
                                           int minInteriorPoints, double pdfPathCm, const char* histTitle = nullptr);

  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke, Float_t mass);
  DEdxFreeRangeFitResult dEdxLikelihoodFreeRangeFit(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0,
                                                    double Lmax, double step, double measuredTrackLengthCm,
                                                    bool computeMomentum, double pdfPathCm = 0.65,
                                                    bool computeDedxBiasDiagnostics = true,
                                                    double scanStepFineCm = 0.,
                                                    double lowPMomentumRefineGeV = 0.2);

  /// For π (PDG 211): same free-range scan as GetdEdxLikelihoodFreeRangeFit, map each trial offset L to
  /// momentum via R_eff = measured_length + L and RangeCmToMomentumGeV. Duplicate p bins keep max log L.
  /// Returns false if fewer than two usable (p, logL) points.
  /// Same interior selection as π free-range TLE: collection hits, optional RR cap, finite dE/dx and RR, RR > 0.
  bool InteriorPionCollectionPlaneDedxRr(AnaParticlePD* part, double maxResidualRangeCm, int minInteriorPoints,
                                          int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                          double dedxMaxMeVcm, std::vector<double>& dedxMeVcmOut,
                                          std::vector<double>& rrCmOut);

  /// R_eff baseline for free-range momentum mapping: particle Length when valid, else max interior RR.
  double MeasuredTrackLengthCmForPionFreeRange(const AnaParticlePD* part, const std::vector<double>& rrInteriorCm);

  bool BuildPionFreeRangeLogLikelihoodVsMomentumCurve(AnaParticlePD* part, double Lmax, double step,
                                                       int minInteriorPoints, int skipHitsFirst, int skipHitsLast,
                                                       double dedxMinMeVcm, double dedxMaxMeVcm, double pdfPathCm,
                                                       std::vector<double>& pGeV, std::vector<double>& logL,
                                                       double scanStepFineCm = 0.,
                                                       double lowPMomentumRefineGeV = 0.2);
  bool BuildPionFreeRangeLogLikelihoodVsMomentumCurveFromVectors(
      const std::vector<double>& rrCm, const std::vector<double>& dedxMeVcm, double measuredTrackLengthCm, double Lmax,
      double step, int minInteriorPoints, int skipFirst, int skipLast, double dedxMinMeVcm, double dedxMaxMeVcm,
      double pdfPathCm, std::vector<double>& pGeV, std::vector<double>& logL, double scanStepFineCm = 0.,
      double lowPMomentumRefineGeV = 0.2);

  /// Pion (PDG 211): mean -log PDF per collection hit in the Bragg window (0 < RR <= maxBraggResidualRangeCm), using
  /// measured RR with no free-range length shift (L = 0). Same Landau/Vavilov PDF as TLE free-range fit.
  /// Smaller meanNegLogL indicates better agreement with the stopping-pion dE/dx template in that window.
  /// Returns false if fewer than minBraggHits hits contribute; meanNegLogL is set to -999 on failure.
  bool ComputePionBraggWindowMeanNegLogLikelihoodVsTemplate(AnaParticlePD* part, double maxBraggResidualRangeCm,
                                                            int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                                            double dedxMaxMeVcm, double pdfPathCm, int minBraggHits,
                                                            double& meanNegLogL, int& nHitsUsed);

  /// Pion (PDG 211): χ²_π± mean per arXiv:2409.18288 Eq. 6.1 — average over collection hits with 0 < RR < maxResidualRangeCm
  /// of [(dEdx_meas − ⟨dE/dx⟩_BB from Eq. 2.1 at KE(RR))/σ]²; KE from measured residual range (CSDA pion curve).
  /// Smaller indicates better agreement with the stopping-pion hypothesis in that window.
  /// If dedxMinMeVcm > 0 and dedxMaxMeVcm > dedxMinMeVcm, apply the same dE/dx bounds as optional robustness (paper §6.1
  /// does not require them). skipHitsFirst/Last trim hit indices before the RR window (default 0 = paper-literal).
  bool ComputePionBraggWindowChi2PiEq61(AnaParticlePD* part, double maxResidualRangeCm, double sigmaDedxMeVcm,
                                        int minHits, int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                        double dedxMaxMeVcm, double& meanChi2, int& nHitsUsed);

  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
  double KineticEnergyMeVFromResidualRangeCm(TGraph* tg_ke, double range_cm);
  double GetWmax(double KE, double mass);
  double GetLandauXi(double KE, double dx, double mass);
  double dEdxPDF(double *x, double *par);

}

#endif
