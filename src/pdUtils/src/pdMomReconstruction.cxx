#include "pdMomReconstruction.hxx"
#include <TVector3.h>
#include <iostream>
#include <cmath>

// Calorimetric momentum reconstruction for interacting pions
// Sums total deposited energy from particle + all daughters

namespace pdMomReconstruction {

//***************************************************************
double CalculateDepositedEnergy(AnaParticlePD* particle, int plane) {
//***************************************************************

  if (!particle) {
    return -999.0;
  }

  double totalDeposited = 0.0;
  int nValidHits = 0;

  // Sum energy from all hits on this plane
  for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
    AnaHitPD& hit = particle->Hits[plane][i];

    // Skip invalid hits
    if (hit.dEdx <= 0 || hit.dEdx > 1000 || hit.dEdx == -999) {
      continue;
    }

    // Use pitch for accurate path length through detector
    double pathLength = 0.0;
    if (hit.Pitch > 0 && hit.Pitch != -999) {
      pathLength = hit.Pitch;
    } else {
      // Fallback: use nominal pitch for collection plane
      pathLength = 0.4792; // cm
    }

    // Energy deposited at this hit: dE = dE/dx * dx
    // dE/dx is in MeV/cm, pathLength is in cm → dE in MeV
    totalDeposited += hit.dEdx * pathLength;
    nValidHits++;
  }

  // Require minimum number of hits for reliability
  if (nValidHits < 3) {
    return -999.0;
  }

  return totalDeposited; // MeV
}

//***************************************************************
double EnergyToMomentum(double totalEnergyMeV, double massMeV) {
//***************************************************************

  // Relativistic energy-momentum relation:
  // E² = (pc)² + (mc²)²
  // p = sqrt(E² - m²) / c

  if (totalEnergyMeV <= massMeV) {
    return -999.0; // Unphysical: total energy must be at least rest mass
  }

  double E2 = totalEnergyMeV * totalEnergyMeV;
  double m2 = massMeV * massMeV;
  double p2 = E2 - m2;

  if (p2 < 0) {
    return -999.0;
  }

  double momentum_MeV = sqrt(p2);
  double momentum_GeV = momentum_MeV / 1000.0; // Convert to GeV/c

  return momentum_GeV;
}

//***************************************************************
double GetRestMass(int pdg) {
//***************************************************************

  // Return rest mass in MeV for common particles
  int absPdg = abs(pdg);

  switch(absPdg) {
    case 211:  // Pion (charged)
      return 139.57;
    case 111:  // Pion (neutral)
      return 134.98;
    case 2212: // Proton
      return 938.27;
    case 2112: // Neutron
      return 939.57;
    case 11:   // Electron
      return 0.511;
    case 13:   // Muon
      return 105.66;
    case 321:  // Kaon (charged)
      return 493.68;
    case 310:  // K0-short
    case 130:  // K0-long
      return 497.61;
    case 22:   // Gamma
      return 0.0;
    default:
      // Unknown particle - return 0 (conservative)
      return 0.0;
  }
}

//***************************************************************
void CollectAllDescendants(AnaParticlePD* particle, std::vector<AnaParticlePD*>& allDescendants) {
//***************************************************************

  if (!particle) return;

  // Add this particle to the list
  allDescendants.push_back(particle);

  // Recursively add all daughters
  for (size_t i = 0; i < particle->Daughters.size(); i++) {
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(particle->Daughters[i]);
    if (daughter) {
      // Recursive call to collect daughter and all its descendants
      CollectAllDescendants(daughter, allDescendants);
    }
  }
}

//***************************************************************
CalorimetricResult CalculateTotalEnergy(AnaParticlePD* particle,
                                        double particleMass) {
//***************************************************************

  CalorimetricResult result;
  result.totalEnergy = -999;
  result.momentum = -999;
  result.nHitsUsed = 0;
  result.nDaughtersIncluded = 0;
  result.nDescendantsRecursive = 0;
  result.nTerminalParticles = 0;
  result.restMassAdded = 0.0;
  result.valid = false;

  if (!particle) {
    return result;
  }

  int plane = 2; // Collection plane

  // Recursively collect all descendants (particle + children + grandchildren + ...)
  std::vector<AnaParticlePD*> allDescendants;
  CollectAllDescendants(particle, allDescendants);
  result.nDescendantsRecursive = allDescendants.size() - 1; // Don't count the parent itself

  double totalDepositedEnergy = 0.0;
  int totalHits = 0;
  double totalRestMass = 0.0;
  int nTerminal = 0;

  // Process each particle in the descendant tree
  for (size_t i = 0; i < allDescendants.size(); i++) {
    AnaParticlePD* desc = allDescendants[i];
    if (!desc) continue;

    // Calculate deposited energy from this particle
    double descEnergy = CalculateDepositedEnergy(desc, plane);
    bool hasDepositedEnergy = (descEnergy > 0 && descEnergy != -999);

    if (hasDepositedEnergy) {
      totalDepositedEnergy += descEnergy;

      // Count hits
      for (size_t j = 0; j < desc->Hits[plane].size(); j++) {
        if (desc->Hits[plane][j].dEdx > 0 &&
            desc->Hits[plane][j].dEdx != -999 &&
            desc->Hits[plane][j].dEdx < 1000) {
          totalHits++;
        }
      }
    }

    // Check if this is a terminal particle (no descendants or no energy deposition)
    bool hasDescendants = (desc->Daughters.size() > 0);
    bool isTerminal = !hasDescendants || !hasDepositedEnergy;

    if (isTerminal) {
      // Add rest mass of terminal particle using true PDG
      AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(desc->TrueObject);
      if (trueParticle) {
        int truePDG = trueParticle->PDG;
        double restMass = GetRestMass(truePDG);
        totalRestMass += restMass;
        nTerminal++;
      }
    }
  }

  // Count first-level daughters for backward compatibility
  result.nDaughtersIncluded = particle->Daughters.size();

  result.totalEnergy = totalDepositedEnergy;
  result.nHitsUsed = totalHits;
  result.restMassAdded = totalRestMass;
  result.nTerminalParticles = nTerminal;

  // Convert deposited energy to total energy
  // Deposited energy = kinetic energy from all visible particles
  // Add rest masses of terminal descendants (which already includes parent if terminal)
  // Total energy = KE_visible + M_terminal_particles
  // Note: Do NOT add particleMass separately - parent is in allDescendants and counted in totalRestMass if terminal
  double totalEnergy_with_mass = totalDepositedEnergy + totalRestMass;

  // Calculate momentum from total energy
  // p = sqrt(E_total² - M_pion²) where E_total is the reconstructed pion total energy
  result.momentum = EnergyToMomentum(totalEnergy_with_mass, particleMass);

  // Validate result
  result.valid = (result.momentum > 0 &&
                  result.momentum < 10.0 &&    // Reasonable upper limit (GeV/c)
                  result.nHitsUsed >= 5);       // Minimum hits for reliability

  return result;
}

//***************************************************************
Float_t EstimateMomentumCalorimetric(AnaParticlePD* particle, int pdg) {
//***************************************************************

  if (!particle || pdg != 211) {
    // Only implemented for pions currently
    return -999.0;
  }

  // Pion mass in MeV
  const double pionMass = 139.57;

  // Calculate total energy and momentum
  CalorimetricResult result = CalculateTotalEnergy(particle, pionMass);

  // Return momentum if calculation was successful
  if (result.valid) {
    return static_cast<Float_t>(result.momentum);
  }

  return -999.0;
}

} // namespace pdMomReconstruction

