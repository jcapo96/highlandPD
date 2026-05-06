#ifndef pdMomReconstructionMCS_h
#define pdMomReconstructionMCS_h

#include "pdDataClasses.hxx"
#include <TVector3.h>
#include <vector>

namespace pdMomReconstruction {

/// Segmentation and Highland likelihood knobs for pion MCS (reco or joint axis).
struct MCSLikelihoodConfig {
  double radiationLengthCm = 14.0;
  double targetSegmentLengthCm = 10.0;
  double minSegmentLengthCm = 0.5;
  double theta0FloorRad = 1e-6;
  double maxAbsDeltaThetaRad = -1.0;
  bool useDetectorSigma = false;
  double detectorSigmaFloorRad = 1e-6;
  double detectorSigmaA = 0.0;
  double detectorSigmaC = 0.0;
};

class MCSLikelihood {
public:
  explicit MCSLikelihood(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg = MCSLikelihoodConfig());

  double ComputeNLL(double momentumGeV) const;
  bool HasSamples() const { return !theta_xz_.empty(); }
  size_t SampleCount() const { return theta_xz_.size(); }

private:
  std::vector<double> theta_xz_;
  std::vector<double> theta_yz_;
  std::vector<double> x_over_x0_;
  std::vector<double> cumulative_length_cm_;
  double theta0_floor_rad_ = 1e-6;
  double radiation_length_cm_ = 14.0;
  bool use_detector_sigma_ = false;
  double detector_sigma_floor_rad_ = 1e-6;
  double detector_sigma_a_ = 0.0;
  double detector_sigma_c_ = 0.0;
};

struct MCSSegmentGeometry {
  TVector3 startPoint;
  TVector3 endPoint;
  TVector3 centroid;
  TVector3 fittedDirection;
  double arcLengthCm = 0.0;
};

bool BuildPionMcsSegmentGeometry(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                 std::vector<MCSSegmentGeometry>& outSegments);

bool BuildPionMcsSegmentGeometryFromOrderedPositions(const std::vector<TVector3>& orderedPositions,
                                                     const MCSLikelihoodConfig& cfg,
                                                     std::vector<MCSSegmentGeometry>& outSegments);

bool BuildPionMcsScatteringSamples(const AnaParticlePD& track, const MCSLikelihoodConfig& cfg,
                                   std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                   std::vector<double>* rrMidCm = nullptr, std::vector<double>* thetaXZ = nullptr,
                                   std::vector<double>* thetaYZ = nullptr,
                                   std::vector<double>* segmentFitChi2Ndf = nullptr,
                                   std::vector<double>* segment1RrToEndCm = nullptr,
                                   std::vector<double>* segment2RrToEndCm = nullptr,
                                   std::vector<double>* segmentLengthCm = nullptr,
                                   std::vector<double>* segmentFitChi2NdfSingle = nullptr,
                                   std::vector<double>* segmentRrToEndCmSingle = nullptr);

bool BuildPionMcsScatteringSamplesFromOrderedPositions(const std::vector<TVector3>& orderedPositions,
                                                       const MCSLikelihoodConfig& cfg,
                                                       std::vector<double>& deltaTheta, std::vector<double>& xOverX0,
                                                       std::vector<double>* rrMidCm = nullptr,
                                                       std::vector<double>* thetaXZ = nullptr,
                                                       std::vector<double>* thetaYZ = nullptr,
                                                       std::vector<double>* segmentFitChi2Ndf = nullptr,
                                                       std::vector<double>* segment1RrToEndCm = nullptr,
                                                       std::vector<double>* segment2RrToEndCm = nullptr,
                                                       std::vector<double>* segmentLengthCm = nullptr,
                                                       std::vector<double>* segmentFitChi2NdfSingle = nullptr,
                                                       std::vector<double>* segmentRrToEndCmSingle = nullptr);

bool BuildPionMCSLogLikelihoodVsMomentumCurve(const AnaParticlePD& track, const std::vector<double>& pAxisGeV,
                                              const MCSLikelihoodConfig& cfg, std::vector<double>& logL);

double ComputePionMcsNLLFromScatteringSamples(const std::vector<double>& theta_xz, const std::vector<double>& theta_yz,
                                                const std::vector<double>& x_over_x0,
                                                const std::vector<double>& cumulative_length_cm_midtrack,
                                                double momentumGeV, const MCSLikelihoodConfig& cfg);

/// MCS likelihood / segmentation settings plus scan and segment-drop controls (pionMomentumAnalysis MCS parameters).
struct PionMCSConfig {
  MCSLikelihoodConfig likelihood{};
  int dropFirstNSegments = 0;
  int dropLastNSegments = 0;
  double momentumScanMaxGeV = 2.5;
};

struct PionMCSResult {
  double bestMomentumGeV = -999.0;
  double bestLogLikelihood = -999.0;
  bool valid = false;
};

struct PionMCSRecoBuffers {
  void clear() {
    deltaThetaPair.clear();
    deltaThetaXz.clear();
    deltaThetaYz.clear();
    segmentLengthPair.clear();
    segmentLengthSingle.clear();
    segmentFitChi2NdfSingle.clear();
    segmentRrToEndSingle.clear();
    scatterUsedPair.clear();
    scatterUsedXz.clear();
    scatterUsedYz.clear();
    scatterUsedRR.clear();
  }
  std::vector<double> deltaThetaPair;
  std::vector<double> deltaThetaXz;
  std::vector<double> deltaThetaYz;
  std::vector<double> segmentLengthPair;
  std::vector<double> segmentLengthSingle;
  std::vector<double> segmentFitChi2NdfSingle;
  std::vector<double> segmentRrToEndSingle;
  std::vector<double> scatterUsedPair;
  std::vector<double> scatterUsedXz;
  std::vector<double> scatterUsedYz;
  std::vector<double> scatterUsedRR;
};

struct PionMCSTrueBuffers {
  void clear() {
    deltaThetaPair.clear();
    deltaThetaXz.clear();
    deltaThetaYz.clear();
    segmentLengthPair.clear();
    segmentLengthSingle.clear();
    segmentFitChi2NdfSingle.clear();
    segmentRrToEndSingle.clear();
  }
  std::vector<double> deltaThetaPair;
  std::vector<double> deltaThetaXz;
  std::vector<double> deltaThetaYz;
  std::vector<double> segmentLengthPair;
  std::vector<double> segmentLengthSingle;
  std::vector<double> segmentFitChi2NdfSingle;
  std::vector<double> segmentRrToEndSingle;
};

bool EstimatePionMomentumMCSReco(AnaParticlePD* track, const PionMCSConfig& cfg, PionMCSResult& out,
                                 PionMCSRecoBuffers& reco, std::vector<double>* pAxisGeVOut = nullptr,
                                 std::vector<double>* logLOut = nullptr);

bool EstimatePionMomentumMCSTrue(const std::vector<TVector3>& orderedScePositions, const PionMCSConfig& cfg,
                                 PionMCSResult& out, PionMCSTrueBuffers& buf);

}  // namespace pdMomReconstruction

#endif
