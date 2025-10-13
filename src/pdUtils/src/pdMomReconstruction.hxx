#ifndef pdMomReconstruction_h
#define pdMomReconstruction_h

// Calorimetric momentum reconstruction for interacting pions
// Sums total deposited energy from particle + all daughters using pitch-corrected path lengths

#include "pdDataClasses.hxx"

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

  // Main function: estimate momentum from deposited energy (particle + daughters)
  Float_t EstimateMomentumCalorimetric(AnaParticlePD* particle, int pdg = 211);

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

#endif

