#ifndef pdAnalysisUtils_h
#define pdAnalysisUtils_h

#include "ParticleId.hxx"
#include "AnalysisUtils.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "pdDataClasses.hxx"
#include "TProfile.h"


namespace pdAnaUtils{

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

  // Create reconstructed vertex candidates from particles
  std::vector<AnaVertexPD*> CreateReconstructedVertices(AnaEventB& event, double maxVertexRadius = 30.0, double maxDaughterDistance = 5.0);

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

  Float_t GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG);
  Float_t GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR);
  Float_t dEdxLikelihood(TGraph* tg, TGraph* tg_ke,
			 Float_t mass);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR);
  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke,
				  Float_t mass);
  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
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
  /// The function fits a 3D line to hits within trackLength cm of the position returned by DefinePosition
  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition);

  /// Overloaded version for backward compatibility (defaults to start position)
  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength = 15.0);

  /// Extrapolate track for true particles
  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition);

  /// Overloaded version of ExtrapolateTrack for true particles (defaults to start position)
  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength = 15.0);

  /// Helper function to fit a line to a set of 3D points using PCA
  /// @param points Vector of 3D points to fit
  /// @param fitParams Output: Vector with 6 parameters [x0, y0, z0, dx, dy, dz] (point + direction)
  void FitLineToPoints(const std::vector<TVector3>& points, std::vector<double>& fitParams);
  void FitLineToPointsPCA(const std::vector<TVector3>& points, std::vector<double>& fitParams);

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

  /// Create neutral particles from vertices by checking particles within sphere and impact parameter
  /// @param event The event containing particles
  /// @param vertices Vector of vertices to check against
  /// @param VertexRadius Radius of sphere around vertex to check for particle end positions
  /// @param ImpactParameter Maximum impact parameter for neutral particle creation
  /// @return Vector of created AnaNeutralParticlePD objects
  std::vector<AnaNeutralParticlePD*> CreateAnaNeutralParticles(AnaEventB& event, const std::vector<AnaVertexPD*>& vertices, double VertexRadius, double ImpactParameter);

  /// Find vertex position by fitting lines to daughter particles and finding closest points
  /// @param vertex The vertex to find position for
  void FindVertexPosition(AnaVertexPD* vertex);

  /// Find vertex position by fitting lines to true daughter particles and finding closest points
  /// @param vertex The true vertex to find position for
  void FindVertexPosition(AnaTrueEquivalentVertexPD* vertex);
}

#endif
