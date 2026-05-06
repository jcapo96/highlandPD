#ifndef pdJointK0sPionMomentum_hxx
#define pdJointK0sPionMomentum_hxx

#include "pdDataClasses.hxx"
#include "TVector3.h"
#include <limits>
#include <vector>

class TH2F;
struct JointK0sPionMomentumGridResult {
  bool ok = false;
  Float_t p1 = -999.f;
  Float_t p2 = -999.f;
  Float_t bestScore = -999.f;
  Float_t invMassAtBest = -999.f;
  Float_t constraintRatioR = -999.f;
};

struct JointK0sSigmaMEventDiagnostics {
  double sigma_p1_gev = std::numeric_limits<double>::quiet_NaN();
  double sigma_p2_gev = std::numeric_limits<double>::quiet_NaN();
  double sigma_m_event_gev = std::numeric_limits<double>::quiet_NaN();
  double dm_dp1 = std::numeric_limits<double>::quiet_NaN();
  double dm_dp2 = std::numeric_limits<double>::quiet_NaN();
  double cos_theta = std::numeric_limits<double>::quiet_NaN();
};

namespace pdJointK0sPionMomentum {

double InterpolateLogLikelihoodClamped(const std::vector<double>& pAxis, const std::vector<double>& logL, double p);

bool PionPairInvariantMassGeV(double a1, double a2, const TVector3& dir1, const TVector3& dir2, double& massGeV);

bool PionPiPiInvariantMassAndDerivativesGeV(double p1, double p2, double cos_theta, double mpi_gev, double& m_gev,
                                           double& dm_dp1, double& dm_dp2);

double EstimateSigmaPFromLogLikelihoodCurve(const std::vector<double>& pAxis, const std::vector<double>& logL,
                                            double p0_gev);

bool ComputeSigmaMEventGeV(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                           const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                           const TVector3& dir1, const TVector3& dir2, double fallback_sigma_m_gev,
                           JointK0sSigmaMEventDiagnostics& out, double sigmaMinGeV, double sigmaMaxGeV);

JointK0sPionMomentumGridResult FitJointMomentaOnGrid(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                                                     const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                                                     const TVector3& dir1, const TVector3& dir2, double pMinGeV,
                                                     double pMaxGeV, double pStepGeV, double mK0sGeV, double sigmaMassGeV,
                                                     double penaltyScale, int refineFactor);

bool MakeJointK0sObjectiveTH2CoarsePass(const std::vector<double>& p1Axis, const std::vector<double>& logL1,
                                        const std::vector<double>& p2Axis, const std::vector<double>& logL2,
                                        const TVector3& dir1, const TVector3& dir2, double pMinGeV, double pMaxGeV,
                                        double pStepGeV, double mK0sGeV, double sigmaMassGeV, double penaltyScale,
                                        const char* nameObjective, const char* titleObjective, const char* namePenalty,
                                        const char* titlePenalty, const char* nameTrackLL, const char* titleTrackLL,
                                        TH2F*& hObjective, TH2F*& hMassPenalty, TH2F*& hTrackLogLSum);

}  // namespace pdJointK0sPionMomentum

#endif
