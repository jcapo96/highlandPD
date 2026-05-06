#include "pdMomReconstructionFromParams.hxx"
#include "Parameters.hxx"

namespace pdMomReconstruction {

void FillPionTLEFitConfig_FromPionMomentumParams(PionTLEFitConfig& cfg) {
  cfg.minInteriorHits = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxMinInteriorHits");
  cfg.skipHitsFirst = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxSkipHitsFirst");
  cfg.skipHitsLast = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxSkipHitsLast");
  cfg.dedxMinMeVcm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxMinMeVcm");
  cfg.dedxMaxMeVcm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxMaxMeVcm");
  cfg.scanLmaxCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanLmaxCm");
  cfg.scanStepCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanStepCm");
  cfg.scanStepFineCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanStepFineCm");
  cfg.lowPMomentumRefineGeV = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeLowPMomentumRefineGeV");
  cfg.dedxPdfPathCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxPdfPathCm");
}

void FillPionMCSConfig_FromPionMomentumParams(PionMCSConfig& cfg) {
  cfg.likelihood.radiationLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSRadiationLengthCm");
  cfg.likelihood.targetSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSTargetSegmentLengthCm");
  cfg.likelihood.minSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSMinSegmentLengthCm");
  cfg.likelihood.theta0FloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCStheta0FloorRad");
  cfg.likelihood.maxAbsDeltaThetaRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSMaxAbsDeltaThetaRad");
  cfg.likelihood.useDetectorSigma = ND::params().GetParameterI("pionMomentumAnalysis.MCSUseDetectorSigma") != 0;
  cfg.likelihood.detectorSigmaFloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaFloorRad");
  cfg.likelihood.detectorSigmaA = ND::params().HasParameter("pionMomentumAnalysis.MCSDetectorSigmaA")
                                      ? ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaA")
                                      : 0.0;
  cfg.likelihood.detectorSigmaC = ND::params().HasParameter("pionMomentumAnalysis.MCSDetectorSigmaC")
                                      ? ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaC")
                                      : 0.0;
  cfg.dropFirstNSegments = ND::params().GetParameterI("pionMomentumAnalysis.MCSDropFirstNSegments");
  cfg.dropLastNSegments = ND::params().GetParameterI("pionMomentumAnalysis.MCSDropLastNSegments");
  cfg.momentumScanMaxGeV = ND::params().HasParameter("pionMomentumAnalysis.MCSMomentumScanMaxGeV")
                               ? ND::params().GetParameterD("pionMomentumAnalysis.MCSMomentumScanMaxGeV")
                               : 2.5;
}

void FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams(JointK0sTwoPionGridFitConfig& cfg) {
  cfg.tle.minInteriorHits =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
          ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
          : 15;
  cfg.tle.skipHitsFirst = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
                              ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
                              : 3;
  cfg.tle.skipHitsLast = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
                             ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
                             : 3;
  cfg.tle.dedxMinMeVcm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
                           ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
                           : 0.5;
  cfg.tle.dedxMaxMeVcm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
                           ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
                           : 5.0;
  cfg.tle.scanLmaxCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                           ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                           : 450.;
  cfg.tle.scanStepCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepCm")
                           ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepCm")
                           : 1.;
  cfg.tle.scanStepFineCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepFineCm")
                               ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepFineCm")
                               : 0.;
  cfg.tle.lowPMomentumRefineGeV =
      ND::params().HasParameter("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
          ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeLowPMomentumRefineGeV")
          : 0.2;
  cfg.tle.dedxPdfPathCm = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                              ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
                              : 0.65;

  cfg.mcs.likelihood.radiationLengthCm =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSRadiationLengthCm")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSRadiationLengthCm")
          : 14.0;
  cfg.mcs.likelihood.minSegmentLengthCm =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSMinSegmentCm")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSMinSegmentCm")
          : 0.5;
  cfg.mcs.likelihood.theta0FloorRad =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSTheta0FloorRad")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSTheta0FloorRad")
          : 1e-6;
  cfg.mcs.likelihood.maxAbsDeltaThetaRad =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSMaxAbsDeltaThetaRad")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSMaxAbsDeltaThetaRad")
          : -1.0;
  cfg.mcs.likelihood.useDetectorSigma =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSUseDetectorSigma") &&
      ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMCSUseDetectorSigma") != 0;
  cfg.mcs.likelihood.detectorSigmaFloorRad =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaFloorRad")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaFloorRad")
          : 1e-6;
  cfg.mcs.likelihood.detectorSigmaA =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaA")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaA")
          : 0.0;
  cfg.mcs.likelihood.detectorSigmaC =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSDetectorSigmaC")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSDetectorSigmaC")
          : 0.0;
  cfg.mcs.dropFirstNSegments = 0;
  cfg.mcs.dropLastNSegments = 0;
  cfg.mcs.momentumScanMaxGeV = 2.5;

  cfg.pMinGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                    ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                    : 0.05;
  cfg.pMaxGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                    ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                    : 2.0;
  cfg.pStepGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                     ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                     : 0.05;
  cfg.mK0sMassGeV = 0.497611;
  const double sigmaMassMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                  ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                  : 10.0;
  const double sigmaMassMinMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMinMeV")
                                       ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMinMeV")
                                       : 5.0;
  const double sigmaMassMaxMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMaxMeV")
                                       ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMaxMeV")
                                       : 50.0;
  cfg.sigmaMassGeV = sigmaMassMeV * 1e-3;
  cfg.sigmaMassMinGeV = sigmaMassMinMeV * 1e-3;
  cfg.sigmaMassMaxGeV = sigmaMassMaxMeV * 1e-3;
  cfg.massPenaltyScale = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                             ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                             : 1.0;
  cfg.refineFactor = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumRefineFactor")
                         ? ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMomentumRefineFactor")
                         : 2;
  cfg.useMCS = !ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSEnable") ||
             ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMCSEnable") != 0;
  cfg.mcsWeight = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMCSWeight")
                      ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMCSWeight")
                      : 1.0;
  cfg.useEventSigmaM = !ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") ||
                       ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") != 0;
}

}  // namespace pdMomReconstruction
