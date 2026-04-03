#ifndef pdAnalysisUtils_h
#define pdAnalysisUtils_h

#include "ParticleId.hxx"
#include "AnalysisUtils.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "pdDataClasses.hxx"
#include "TProfile.h"

class TMultiGraph;
class TGraph;

namespace pdAnaUtils{

  /// dE/dx free-range (track-length extension) ML fit: log L, best offset, and optional momentum (GeV/c).
  struct DEdxFreeRangeFitResult {
    Float_t logLikelihood = -999.f;
    Float_t bestOffsetCm = -999.f;
    Float_t momentumGeV = -999.f;
  };

  /// Computes the range momentum
  Float_t ComputeRangeMomentum(double trkrange, int pdg);

  /// Computes the CSDARange
  Float_t ComputeCSDARange(double beammom, int pdg);

  /// Computes the kinetic energy
  Float_t ComputeKineticEnergy(const AnaParticlePD &part);

  // Compute distances between daughters and vertex
  void ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance);

  /// Find the beam true particle
  AnaTrueParticle* GetBeamTrueParticle(const AnaSpillB& spill);

  // Get te AnaTrueParticle with a given ID
  AnaTrueParticlePD* GetTrueParticle(AnaEventB* event, Int_t ID);
  AnaTrueParticlePD* GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID);

  AnaParticlePD* GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID);

  // retreieve the BI particle
  AnaParticlePD* GetBeamParticle(const AnaEventC& event);

  // retreieve the true BI particle
  AnaTrueParticlePD* GetTrueBeamParticle(const AnaEventC& event);

  //checks if the beam particle selected by Pandora has been correctly selected. Basic implementation for the moment
  bool isBeamLike(AnaParticlePD* part, AnaBeamPD* beam);

  // Compute the PID chi2 and ndf depending of pdg hypothesis and up to a given residual range
  std::pair< double, int > Chi2PID(const AnaParticlePD& part, const int pdg);
  std::pair< double, int > Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR);

  // Compute the track length using the hit positions
  Float_t ComputeTrackLengthFromHitPosition(const AnaParticlePD* part);
  Float_t ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part);

  void ComputeParticlePositionAndDirection(AnaParticlePD* part);

  // Compute the truncated mean of an std vector
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> dEdx);
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD> hits);

  // Distance mother daughter
  Float_t ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);

  // Cos mother daughter
  Float_t ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);

  // Average dEdx
  Float_t ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange = 9999);

  bool IsStoppingInFV(AnaParticlePD *part);

  int GetHitTPCid(AnaHitPD& hit);
  int GetPosTPCid(TVector3 pos);

  void EstimateHitsDirection(AnaParticlePD* part);

  void ComputeResidualRange(AnaParticlePD* part);

  Double_t ComputeDepositedEnergy(AnaParticlePD* part);

  Double_t EstimateTrueMomAtAPABorder(AnaParticlePD* part);

  /// Calculate momentum from calorimetric energy, optionally including decay products
  /// @param part The particle to calculate momentum for
  /// @param pdg The PDG code of the particle hypothesis
  /// @param includeDecayProducts If true, add energy from daughter particles
  /// @return Momentum in MeV/c
  Float_t ComputeCalorimetricMomentum(AnaParticlePD* part, int pdg, bool includeDecayProducts = false);

  // Compute invariant mass for true particles (helper function)
  Float_t ComputeTrueInvariantMass(const AnaTrueParticlePD& part1, const AnaTrueParticlePD& part2, Float_t mass1, Float_t mass2);

  Double_t ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts);

  void GetBeamQualityCuts(AnaEventPD* event,
			  double &mean_x, double &mean_y, double &mean_z,
			  double &sigma_x, double &sigma_y, double &sigma_z,
			  double &cos);

  Float_t GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG,
                            double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  Float_t GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                   double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  Float_t dEdxLikelihood(TGraph* tg, TGraph* tg_ke,
			 Float_t mass);
  /// @param landauTruncMinRRCm Only hits with residual range >= this (cm) may be removed as tail hits; low-RR hits are kept.
  /// @param landauTailHitDropFraction Fraction of those high-RR hits (by largest dE/dx) to drop; 0 disables.
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG,
                                                      double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                              double landauTruncMinRRCm = 0.,
                                                              double landauTailHitDropFraction = 0.);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax = 500.,
                                                     double step = 0.5, double landauTruncMinRRCm = 0.,
                                                     double landauTailHitDropFraction = 0.);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                              double Lmax = 500., double step = 0.5,
                                                              double landauTruncMinRRCm = 0.,
                                                              double landauTailHitDropFraction = 0.);

  /// Overlay for diagnostics: measured dE/dx vs RR, same dE/dx vs (RR+L_best), and mean dE/dx vs RR from
  /// dedx_range_pi template. Same interior-hit sample and L-scan as GetdEdxLikelihoodFreeRangeFit (does not affect
  /// momentum assigned elsewhere beyond shared algorithms if you call both).
  /// Caller owns the returned TMultiGraph and should delete it or write it to a TFile.
  /// @param xAxisTitle Optional full X-axis label; if null, uses "Residual range [cm]".
  TMultiGraph* MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax = 500., double step = 0.5,
                                                  double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.,
                                                  const char* xAxisTitle = nullptr);

  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke,
				  Float_t mass);
  DEdxFreeRangeFitResult dEdxLikelihoodFreeRangeFit(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0, double Lmax,
                                                   double step, double measuredTrackLengthCm, bool computeMomentum);
  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
  /// KE [MeV] from ke_vs_range graph (x = residual range [cm]); extrapolates when ROOT::Eval is non-physical near RR -> 0.
  double KineticEnergyMeVFromResidualRangeCm(TGraph* tg_ke, double range_cm);
  double GetWmax(double KE, double mass);
  double GetLandauXi(double KE, double dx, double mass);
  double dEdxPDF(double *x, double *par);

  //  Float_t ComputeVertexPolarity(AnaParticlePD* part1, AnaParticlePD* part2);

  /// Fit a line to 3D hits within 15 cm of the position defined by DefinePosition
  /// @param part The particle to analyze
  /// @param part Input: Particle to extrapolate
  /// @param fitParams Output: Vector with 6 parameters [x0, y0, z0, dx, dy, dz] for line fit (point + direction)
  /// @param trackLength Input: Length of track used for line fitting (cm)
  /// @param useStartPosition If true, use start position as reference; if false, use end position
  /// @param trackFitDistanceFromStart Arc-length distance from reference where fit window begins (cm)
  /// The function fits a 3D line to hits in [trackFitDistanceFromStart, trackFitDistanceFromStart + trackLength]
  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition,
                        double trackFitDistanceFromStart = 0.0);

  /// Overloaded version for backward compatibility (defaults to start position)
  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength = 15.0);

  /// Extrapolate track for true particles
  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition);

  /// Overloaded version of ExtrapolateTrack for true particles (defaults to start position)
  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength = 15.0);

  /// Define the position to use for calculations (distance, line fitting, etc.)
  /// @param particle The particle to get position from
  /// @param useStartPosition If true, use start position; if false, use end position
  /// @return TVector3 with the position to use for calculations
  /// This function can be modified to return different positions (e.g., extrapolated positions)
  TVector3 DefinePosition(AnaParticlePD* particle, bool useStartPosition);

  /// Overloaded version for backward compatibility (defaults to start position)
  TVector3 DefinePosition(AnaParticlePD* particle);

  /// Find the closest points between two 3D lines
  /// @param line1Params Vector with 6 parameters [x0, y0, z0, dx, dy, dz] for first line
  /// @param line2Params Vector with 6 parameters [x0, y0, z0, dx, dy, dz] for second line
  /// @param closestPoint1 Output: Closest point on first line
  /// @param closestPoint2 Output: Closest point on second line
  /// @return Minimum distance between the two lines
  double FindClosestPointsBetweenLines(const std::vector<double>& line1Params,
                                      const std::vector<double>& line2Params,
                                      TVector3& closestPoint1,
                                      TVector3& closestPoint2);

  /// Calculate the impact parameter (distance from a point to a line)
  /// @param lineParams Vector with 6 parameters [x0, y0, z0, dx, dy, dz] for the line
  /// @param point The point to calculate distance from
  /// @return Distance from point to line
  double CalculateImpactParameter(const std::vector<double>& lineParams, const TVector3& point);

}

#endif
