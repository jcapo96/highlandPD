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
}

#endif
