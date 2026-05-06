#ifndef pdMomReconstructionJointK0s_h
#define pdMomReconstructionJointK0s_h

#include "pdDataClasses.hxx"
#include "pdMomReconstructionMCS.hxx"
#include "pdJointK0sPionMomentum.hxx"
#include "pdMomReconstructionTLEFit.hxx"
#include <TVector3.h>
#include <vector>

namespace pdMomReconstruction {

/// Parameters for the K0s two-pion joint grid (same recipe as neutralKaonAnalysis joint block).
struct JointK0sTwoPionGridFitConfig {
  PionTLEFitConfig tle{};
  PionMCSConfig mcs{};
  double pMinGeV = 0.05;
  double pMaxGeV = 2.0;
  double pStepGeV = 0.05;
  double mK0sMassGeV = 0.497611;
  double sigmaMassGeV = 0.01;
  double sigmaMassMinGeV = 0.005;
  double sigmaMassMaxGeV = 0.05;
  double massPenaltyScale = 1.0;
  int refineFactor = 2;
  bool useMCS = true;
  double mcsWeight = 1.0;
  bool useEventSigmaM = true;
};

/// Curves and marginal TLE/MCS modes for post-fit diagnostics (e.g. annihilation χ²); filled when non-null.
struct JointK0sTwoPionFitAuxiliary {
  std::vector<double> p1AxisGeV;
  std::vector<double> logL1Tle;
  std::vector<double> p2AxisGeV;
  std::vector<double> logL2Tle;
  double p1TleArgmaxGeV = -999.;
  double p2TleArgmaxGeV = -999.;
  double p1McsArgmaxGeV = -999.;
  double p2McsArgmaxGeV = -999.;
};

/// TLE curves via \ref EstimatePionMomentumTLEFit, optional MCS on each daughter's TLE momentum axis,
/// event σ_M from marginal TLE curves when enabled, then \ref FitJointMomentaOnGrid.
bool FitJointK0sTwoPionMomentaOnGrid(AnaParticlePD* daughter1, AnaParticlePD* daughter2, const TVector3& dir1,
                                      const TVector3& dir2, const JointK0sTwoPionGridFitConfig& cfg,
                                      JointK0sPionMomentumGridResult& out,
                                      JointK0sSigmaMEventDiagnostics* sigmaDiagOut = nullptr,
                                      JointK0sTwoPionFitAuxiliary* auxOut = nullptr);

/// Shared TLE + MCS log L on extended diagnostic p-axis (neutral kaon multigraphs); matches annihilation joint tuning.
bool BuildNeutralKaonJointDiagnosticsCurvesForDaughter(AnaParticlePD* reco, const PionTLEFitConfig& tle,
                                                      const MCSLikelihoodConfig& mcsLikelihood,
                                                      double diagPMinGeV, std::vector<double>& pAxisDiagOut,
                                                      std::vector<double>& rawLogLTleOut,
                                                      std::vector<double>& rawLogLMcsOut, double& pBestTleOut,
                                                      double& pBestMcsOut);

}  // namespace pdMomReconstruction

#endif
