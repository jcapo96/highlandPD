#ifndef pdMomReconstruction_h
#define pdMomReconstruction_h

// Calorimetric momentum reconstruction for interacting pions
// Sums total deposited energy from particle + all daughters using pitch-corrected path lengths

#include "pdDataClasses.hxx"
#include <TProfile.h>
#include <TGraph.h>
#include <vector>

namespace pdMomReconstruction {

  // Result structure for calorimetric momentum estimation
  struct CalorimetricResult {
    double totalEnergy;          // Total deposited energy (MeV)
    double momentum;             // Calculated momentum (GeV/c)
    int nHitsUsed;              // Total hits from particle + descendants
    int nDaughtersIncluded;     // Number of first-level daughters (deprecated)
    int nDescendantsRecursive;  // Total descendants collected recursively
    int nTerminalParticles;     // Particles at end of chain (no descendants)
    double restMassAdded;       // Total rest mass added from terminal particles (MeV)
    bool valid;                 // Whether calculation succeeded
  };

  // Result structure for track-length extension fit
  struct ExtensionFitResult {
    double extension;        // Track extension length (cm)
    double effectiveRange;   // Total effective range (cm)
    double momentum;         // Estimated momentum (GeV/c)
    double chi2;             // Chi-squared of fit
    int ndf;                 // Number of degrees of freedom
    bool valid;              // Whether fit was successful
  };

  // Result structure for likelihood-based momentum estimation
  struct LikelihoodFitResult {
    double momentum;           // Best-fit momentum (GeV/c)
    double logLikelihood;      // Log-likelihood at best fit
    double uncertainty;        // Estimated uncertainty (GeV/c)
    int nHitsUsed;             // Number of hits used in fit
    bool valid;                // Whether fit was successful
  };

  // Calorimetric method: estimate momentum from deposited energy (particle + daughters)
  Float_t EstimateMomentum(AnaParticlePD* particle, int pdg = 211);

  // Track-length extension method
  Float_t EstimateMomentumByExtension(AnaParticlePD* particle, int pdg = 211);
  ExtensionFitResult EstimateMomentumDetailedByExtension(AnaParticlePD* particle, int pdg = 211);
  Float_t EstimateMomentumFromRange(AnaParticlePD* particle);
  Float_t EstimatePionMomentumFromCSDA(const AnaParticlePD* particle);
  Float_t EstimatePionKineticEnergyFromCSDA(const AnaParticlePD* particle);

  // Likelihood method
  Float_t EstimateMomentumByLikelihood(AnaParticlePD* particle, int pdg = 211);

  // Extension helpers
  ExtensionFitResult FitTrackLengthExtension(AnaParticlePD* particle,
                                             TProfile* dEdxTemplate,
                                             TGraph* rangeEnergyGraph,
                                             double pionMass);
  double CalculateExtensionChi2(const std::vector<double>& measuredDeDx,
                                const std::vector<double>& measuredRR,
                                double extension,
                                TProfile* dEdxTemplate);
  double RangeToMomentum(double effectiveRange, int pdg, TGraph* rangeEnergyGraph, double mass);

  // Likelihood helpers
  LikelihoodFitResult FitMomentumLikelihood(AnaParticlePD* particle,
                                            TProfile* dEdxTemplate,
                                            TGraph* rangeEnergyGraph,
                                            double particleMass);
  double CalculateLogLikelihood(const std::vector<double>& measuredDeDx,
                                const std::vector<double>& distanceFromVertex,
                                const std::vector<double>& pitch,
                                double momentumHypothesis,
                                TProfile* dEdxVsRangeTemplate,
                                TGraph* rangeEnergyGraph,
                                double particleMass);
  double GetExpectedDeDx(double distanceTraveled,
                         double momentum,
                         TProfile* dEdxVsRangeTemplate,
                         TGraph* rangeEnergyGraph,
                         double particleMass);
  double MomentumToKineticEnergy(double momentumGeV, double massMeV);
  double KineticEnergyToRange(double kineticEnergyMeV, TGraph* rangeEnergyGraph);
  double RangeToKineticEnergy(double rangeCm, TGraph* rangeEnergyGraph);

  // Calculate total deposited energy for particle + all daughters
  CalorimetricResult CalculateTotalEnergy(AnaParticlePD* particle, double particleMass);

  // Calculate deposited energy for a single particle using pitch
  double CalculateDepositedEnergy(AnaParticlePD* particle, int plane = 2);

  // Convert total energy to momentum given particle mass
  // E² = p²c² + m²c⁴ → p = sqrt(E² - m²) / c
  double EnergyToMomentum(double totalEnergyMeV, double massMeV);

  // Recursively collect all descendants (children, grandchildren, etc.)
  void CollectAllDescendants(AnaParticlePD* particle, std::vector<AnaParticlePD*>& allDescendants);

  // Get rest mass in MeV for a given PDG code
  double GetRestMass(int pdg);

} // namespace pdMomReconstruction

namespace pdMomShared {

  bool LoadPionTemplates(TProfile*& pionDEdxTemplate, TGraph*& pionRangeEnergyGraph,
                         const char* graphCloneName = "pionRangeEnergyGraph_shared");

  double MomentumToKineticEnergyMeV(double momentumGeV, double massMeV);
  double KineticEnergyMeVToRangeCm(double kineticEnergyMeV, TGraph* rangeEnergyGraph);
  double RangeCmToKineticEnergyMeV(double rangeCm, TGraph* rangeEnergyGraph);
  double RangeCmToMomentumGeV(double rangeCm, int pdg, TGraph* rangeEnergyGraph, double massMeV);

}

#include "pdMomReconstructionTLEFit.hxx"
#include "pdMomReconstructionMCS.hxx"
#include "pdMomReconstructionJointK0s.hxx"
#include "pdMomReconstructionFromParams.hxx"

#endif

