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

  // retreieve the BI particle
  AnaParticlePD* GetBeamParticle(const AnaEventC& event);

  // retreieve the true BI particle
  AnaTrueParticlePD* GetTrueBeamParticle(const AnaEventC& event);

  //checks if the beam particle selected by Pandora has been correctly selected. Basic implementation for the moment
  bool isBeamLike(AnaParticlePD* part, AnaBeamPD* beam);
  
  // Compute the PID chi2 and ndf for protons
  std::pair< double, int > Chi2PID(const AnaParticlePD& part, const int pdg);

  // Compute the track length using the hit positions
  Float_t ComputeTrackLengthFromHitPosition(const AnaParticlePD* part); 
  Float_t ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part); 

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
}

#endif
