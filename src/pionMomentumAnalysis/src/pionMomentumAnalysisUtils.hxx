#ifndef pionMomentumAnalysisUtils_h
#define pionMomentumAnalysisUtils_h

#include "pdDataClasses.hxx"
#include "OutputManager.hxx"

#include <vector>

namespace pionMomentumAnaUtils {

  struct MomentumDiagConfig {
    bool enableMomentumDiagnosticMultigraphs = false;
    bool freeRangeComputeDedxBiasDiagnostics = false;
    int freeRangeDedxMinInteriorHits = 15;
    int freeRangeDedxSkipHitsFirst = 3;
    int freeRangeDedxSkipHitsLast = 3;
    double freeRangeDedxMinMeVcm = 0.5;
    double freeRangeDedxMaxMeVcm = 5.0;
    double freeRangeScanLmaxCm = 450.0;
    double freeRangeScanStepCm = 1.0;
    double freeRangeDedxPdfPathCm = 0.65;
    double mcsRadiationLengthCm = 14.0;
  };

  void AddCustomCategories();
  /// Truth-only stopping vs threshold on true |p| at end: 1 = stopping, 2 = not, 3 = truth but bad end momentum, else CATNOTRUTH.
  int MainTrueStoppingCode(AnaParticlePD* mainTrack, double stoppingMaxTrueEndMomentumGeV);
  int FillMainTrueStoppingCategory(AnaParticlePD* mainTrack, double stoppingMaxTrueEndMomentumGeV);
  void MaybeAccumulateMainTrackMomentumDiagnostics(AnaParticlePD* mainTrack, const AnaEventB& event, int stoppingCode,
                                                   bool runTLE,
                                                   const MomentumDiagConfig& cfg, double bestTleMomentumGeV,
                                                   double bestMcsMomentumGeV, const std::vector<double>& pAxisTleGeV,
                                                   const std::vector<double>& logLTle,
                                                   const std::vector<double>& pAxisMcsGeV,
                                                   const std::vector<double>& logLMcs,
                                                   const std::vector<double>& mcsDeltaTheta,
                                                   const std::vector<double>& mcsSegmentLengthCm,
                                                   const std::vector<double>& mcsDeltaThetaKept,
                                                   const std::vector<double>& mcsResidualRangeKeptCm,
                                                   const std::vector<double>* mcsDeltaThetaTrue = nullptr,
                                                   const std::vector<double>* mcsSegmentLengthTrueCm = nullptr);
  void WriteMomentumDiagnostics(OutputManager& output);

}

#endif
