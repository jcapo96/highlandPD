#include "pdMomReconstruction.hxx"
#include "pdAnalysisUtils.hxx"
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// Calorimetric momentum reconstruction for interacting pions
// Sums total deposited energy from particle + all daughters

namespace pdMomReconstruction {

namespace {
static constexpr double kPionMassMeV = 139.57;
} // namespace

namespace {
// CSDA baseline curve for pions from mean Bethe–Bloch dE/dx:
// X = range [cm], Y = KE [MeV].
//
// Paper reference: E_K^{full} from CSDA using Eq. (2.1).
TGraph* PionCsdaRangeToKeGraphBB() {
  static TGraph* tg = nullptr;
  if (tg) return tg;

  // KE grid (log-spaced) where we tabulate CSDA range.
  constexpr double kKeMinMeV = 1.0;     // avoid beta->0 singularities
  constexpr double kKeMaxMeV = 2000.0; // sufficient for this analysis
  constexpr int kNKeGrid = 80;
  constexpr int kNIntegrationSteps = 500;

  std::vector<std::pair<double, double>> pts; // (range_cm, ke_mev)
  pts.reserve(kNKeGrid);

  for (int i = 0; i < kNKeGrid; ++i) {
    const double t = (kNKeGrid > 1) ? static_cast<double>(i) / (kNKeGrid - 1) : 0.;
    const double keMeV = kKeMinMeV * std::exp(t * std::log(kKeMaxMeV / kKeMinMeV));

    // CSDA: range = ∫ dE / (dE/dx(E)), with dE/dx from Bethe–Bloch.
    // Midpoint rule in energy.
    const double eMin = kKeMinMeV;
    const double eMax = keMeV;
    if (!(eMax > eMin) || !std::isfinite(eMax)) continue;

    const double dE = (eMax - eMin) / static_cast<double>(kNIntegrationSteps);
    if (!(dE > 0.) || !std::isfinite(dE)) continue;

    double rangeCm = 0.;
    for (int j = 0; j < kNIntegrationSteps; ++j) {
      const double eMid = eMin + (static_cast<double>(j) + 0.5) * dE;
      const double dedx = pdAnaUtils::GetdEdxBetheBloch(eMid, kPionMassMeV); // MeV/cm
      if (!(dedx > 0.) || !std::isfinite(dedx)) continue;
      rangeCm += dE / dedx;
    }

    if (rangeCm > 0. && std::isfinite(rangeCm)) {
      pts.emplace_back(rangeCm, keMeV);
    }
  }

  if (pts.empty()) return nullptr;

  std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

  tg = new TGraph(static_cast<int>(pts.size()));
  tg->SetName("pion_csda_ke_vs_range_bb");
  for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
    tg->SetPoint(i, pts[i].first, pts[i].second); // X: range_cm, Y: KE_MeV
  }
  return tg;
}
} // namespace

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
Float_t EstimateMomentum(AnaParticlePD* particle, int pdg) {
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

//***************************************************************
double CalculateExtensionChi2(const std::vector<double>& measuredDeDx,
                              const std::vector<double>& measuredRR,
                              double extension,
                              TProfile* dEdxTemplate) {
//***************************************************************

  if (!dEdxTemplate || measuredDeDx.size() != measuredRR.size() || measuredDeDx.empty()) {
    return 9999.0;
  }

  double chi2 = 0.0;
  int nPoints = 0;

  for (size_t i = 0; i < measuredDeDx.size(); i++) {
    const double effectiveRR = measuredRR[i] + extension;
    const double expectedDeDx = dEdxTemplate->Interpolate(effectiveRR);

    if (expectedDeDx <= 0 || expectedDeDx > 1000) {
      continue;
    }

    double sigma = 0.04231 + 0.0001783 * measuredDeDx[i] * measuredDeDx[i];
    sigma *= measuredDeDx[i];
    if (sigma < 0.1) sigma = 0.1;

    const double residual = measuredDeDx[i] - expectedDeDx;
    chi2 += (residual * residual) / (sigma * sigma);
    nPoints++;
  }

  if (nPoints < 3) {
    return 9999.0;
  }

  return chi2;
}

//***************************************************************
double RangeToMomentum(double effectiveRange, int pdg, TGraph* rangeEnergyGraph, double mass) {
//***************************************************************
  return pdMomShared::RangeCmToMomentumGeV(effectiveRange, pdg, rangeEnergyGraph, mass);
}

//***************************************************************
ExtensionFitResult FitTrackLengthExtension(AnaParticlePD* particle,
                                           TProfile* dEdxTemplate,
                                           TGraph* rangeEnergyGraph,
                                           double pionMass) {
//***************************************************************

  ExtensionFitResult result;
  result.extension = -999;
  result.effectiveRange = -999;
  result.momentum = -999;
  result.chi2 = 9999;
  result.ndf = 0;
  result.valid = false;

  if (!particle || !dEdxTemplate || !rangeEnergyGraph) {
    return result;
  }

  const int plane = 2;
  if (particle->Hits[plane].empty() || particle->Hits[plane].size() < 3) {
    return result;
  }

  std::vector<double> measuredDeDx;
  std::vector<double> distanceFromStart;

  double cumulativeDistance = 0.0;
  TVector3 previousPos = particle->Hits[plane][0].Position;

  for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
    AnaHitPD& hit = particle->Hits[plane][i];

    if (hit.dEdx <= 0 || hit.dEdx > 1000 || hit.dEdx == -999) {
      continue;
    }

    if (i > 0) {
      TVector3 currentPos = hit.Position;
      cumulativeDistance += (currentPos - previousPos).Mag();
      previousPos = currentPos;
    }

    measuredDeDx.push_back(hit.dEdx);
    distanceFromStart.push_back(cumulativeDistance);
  }

  std::vector<double> measuredRR;
  const double totalMeasuredLength = distanceFromStart.empty() ? 0.0 : distanceFromStart.back();
  for (size_t i = 0; i < distanceFromStart.size(); i++) {
    measuredRR.push_back(totalMeasuredLength - distanceFromStart[i]);
  }

  if (measuredDeDx.size() < 3) {
    return result;
  }

  double minChi2 = 9999.0;
  double bestExtension = 0.0;

  const double extensionMin = 0.0;
  const double extensionMax = 500.0;
  const double extensionStep = 0.5;

  for (double ext = extensionMin; ext <= extensionMax; ext += extensionStep) {
    const double chi2 = CalculateExtensionChi2(measuredDeDx, measuredRR, ext, dEdxTemplate);
    if (chi2 < minChi2) {
      minChi2 = chi2;
      bestExtension = ext;
    }
  }

  const double measuredLength = measuredRR.back();
  const double effectiveRange = measuredLength + bestExtension;
  const double momentum = RangeToMomentum(effectiveRange, 211, rangeEnergyGraph, pionMass);

  result.extension = bestExtension;
  result.effectiveRange = effectiveRange;
  result.momentum = momentum;
  result.chi2 = minChi2;
  result.ndf = static_cast<int>(measuredDeDx.size());

  const double chi2ndf = (result.ndf > 0) ? result.chi2 / result.ndf : 9999.0;
  result.valid = (chi2ndf < 5.0 && momentum > 0 && momentum < 10.0);

  return result;
}

//***************************************************************
ExtensionFitResult EstimateMomentumDetailedByExtension(AnaParticlePD* particle, int pdg) {
//***************************************************************
  ExtensionFitResult invalidResult;
  invalidResult.extension = -999;
  invalidResult.effectiveRange = -999;
  invalidResult.momentum = -999;
  invalidResult.chi2 = 9999;
  invalidResult.ndf = 0;
  invalidResult.valid = false;

  if (!particle || pdg != 211) {
    return invalidResult;
  }

  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;
  if (!pdMomShared::LoadPionTemplates(pionDEdxTemplate, pionRangeEnergyGraph,
                                      "pionRangeEnergyGraph_extension")) {
    std::cout << "WARNING: Failed to load pion templates for momentum estimation" << std::endl;
    return invalidResult;
  }

  return FitTrackLengthExtension(particle, pionDEdxTemplate, pionRangeEnergyGraph, kPionMassMeV);
}

//***************************************************************
Float_t EstimateMomentumByExtension(AnaParticlePD* particle, int pdg) {
//***************************************************************
  ExtensionFitResult result = EstimateMomentumDetailedByExtension(particle, pdg);
  if (result.valid) {
    return static_cast<Float_t>(result.momentum);
  }
  return -999.0;
}

//***************************************************************
Float_t EstimateMomentumFromRange(AnaParticlePD* particle) {
//***************************************************************
  if (!particle) return -999.0f;
  if (!(particle->Length > 0.0f) || particle->Length == -999.0f) return -999.0f;

  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;
  if (!pdMomShared::LoadPionTemplates(pionDEdxTemplate, pionRangeEnergyGraph,
                                      "pionRangeEnergyGraph_extension")) {
    return -999.0f;
  }

  const double momentum = RangeToMomentum(static_cast<double>(particle->Length), 211,
                                          pionRangeEnergyGraph, kPionMassMeV);
  if (!(momentum > 0.0) || !std::isfinite(momentum)) return -999.0f;
  return static_cast<Float_t>(momentum);
}

//***************************************************************
Float_t EstimatePionMomentumFromCSDA(const AnaParticlePD* particle) {
//***************************************************************
  if (!particle) return -999.0f;
  if (!(particle->Length > 0.0f) || particle->Length == -999.0f) return -999.0f;

  TGraph* tg = PionCsdaRangeToKeGraphBB();
  if (!tg) return -999.0f;

  const double rangeCm = static_cast<double>(particle->Length);
  const double kineticEnergyMeV = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg, rangeCm);
  if (!(kineticEnergyMeV > 0.0) || !std::isfinite(kineticEnergyMeV)) return -999.f;

  const double momentumMeV =
      std::sqrt(kineticEnergyMeV * kineticEnergyMeV + 2.0 * kPionMassMeV * kineticEnergyMeV);
  return static_cast<Float_t>(momentumMeV / 1000.0);
}

//***************************************************************
Float_t EstimatePionKineticEnergyFromCSDA(const AnaParticlePD* particle) {
//***************************************************************
  if (!particle) return -999.f;
  if (!(particle->Length > 0.0f) || particle->Length == -999.0f) return -999.f;

  TGraph* tg = PionCsdaRangeToKeGraphBB();
  if (!tg) return -999.f;

  const double rangeCm = static_cast<double>(particle->Length);
  const double kineticEnergyMeV = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg, rangeCm);
  if (!(kineticEnergyMeV > 0.0) || !std::isfinite(kineticEnergyMeV)) return -999.f;
  return static_cast<Float_t>(kineticEnergyMeV / 1000.0);
}

//***************************************************************
double MomentumToKineticEnergy(double momentumGeV, double massMeV) {
//***************************************************************
  return pdMomShared::MomentumToKineticEnergyMeV(momentumGeV, massMeV);
}

//***************************************************************
double KineticEnergyToRange(double kineticEnergyMeV, TGraph* rangeEnergyGraph) {
//***************************************************************
  return pdMomShared::KineticEnergyMeVToRangeCm(kineticEnergyMeV, rangeEnergyGraph);
}

//***************************************************************
double RangeToKineticEnergy(double rangeCm, TGraph* rangeEnergyGraph) {
//***************************************************************
  return pdMomShared::RangeCmToKineticEnergyMeV(rangeCm, rangeEnergyGraph);
}

//***************************************************************
double GetExpectedDeDx(double distanceTraveled,
                       double momentumGeV,
                       TProfile* dEdxVsRangeTemplate,
                       TGraph* rangeEnergyGraph,
                       double particleMass) {
//***************************************************************

  if (!dEdxVsRangeTemplate || !rangeEnergyGraph || distanceTraveled < 0) {
    return -999.0;
  }

  double initialKE = MomentumToKineticEnergy(momentumGeV, particleMass);
  if (initialKE <= 0) {
    return -999.0;
  }

  if (distanceTraveled < 1.0) {
    const double totalRange = KineticEnergyToRange(initialKE, rangeEnergyGraph);
    if (totalRange <= 0) return -999.0;
    const double residualRange = totalRange - distanceTraveled;
    if (residualRange < 0) return -999.0;
    return dEdxVsRangeTemplate->Interpolate(residualRange);
  }

  double energyLost = 0.0;
  const double stepSize = 0.5;
  double currentKE = initialKE;

  for (double d = 0; d < distanceTraveled; d += stepSize) {
    currentKE = initialKE - energyLost;
    if (currentKE <= 0) {
      return -999.0;
    }

    const double currentRange = KineticEnergyToRange(currentKE, rangeEnergyGraph);
    if (currentRange <= 0) {
      return -999.0;
    }

    const double currentDeDx = dEdxVsRangeTemplate->Interpolate(currentRange);
    if (currentDeDx <= 0 || currentDeDx > 1000) {
      return -999.0;
    }

    energyLost += currentDeDx * stepSize;
  }

  const double finalKE = initialKE - energyLost;
  if (finalKE <= 0) {
    return -999.0;
  }

  const double finalRange = KineticEnergyToRange(finalKE, rangeEnergyGraph);
  if (finalRange <= 0) {
    return -999.0;
  }

  const double finalDeDx = dEdxVsRangeTemplate->Interpolate(finalRange);
  if (finalDeDx <= 0 || finalDeDx > 1000) {
    return -999.0;
  }

  return finalDeDx;
}

//***************************************************************
double CalculateLogLikelihood(const std::vector<double>& measuredDeDx,
                              const std::vector<double>& distanceFromVertex,
                              const std::vector<double>& pitch,
                              double momentumHypothesis,
                              TProfile* dEdxVsRangeTemplate,
                              TGraph* rangeEnergyGraph,
                              double particleMass) {
//***************************************************************

  if (measuredDeDx.size() != distanceFromVertex.size() ||
      measuredDeDx.size() != pitch.size() ||
      measuredDeDx.empty()) {
    return -9999.0;
  }

  double logLikelihood = 0.0;
  int nValidPoints = 0;
  double totalWeight = 0.0;
  const double nominalPitch = 0.4792;

  for (size_t i = 0; i < measuredDeDx.size(); i++) {
    const double measured = measuredDeDx[i];
    const double distance = distanceFromVertex[i];
    const double hitPitch = pitch[i];

    const double expected = GetExpectedDeDx(distance, momentumHypothesis,
                                            dEdxVsRangeTemplate, rangeEnergyGraph,
                                            particleMass);

    if (expected <= 0 || expected == -999.0) {
      continue;
    }

    double sigma = 0.04231 + 0.0001783 * measured * measured;
    sigma *= measured;
    if (sigma < 0.1) sigma = 0.1;

    double pitchWeight = 1.0;
    if (hitPitch > 0 && hitPitch != -999) {
      const double pitchDeviation = std::fabs(hitPitch - nominalPitch) / nominalPitch;
      pitchWeight = 1.0 / (1.0 + pitchDeviation);
    }

    const double residual = measured - expected;
    const double contribution = -0.5 * (residual * residual) / (sigma * sigma);
    logLikelihood += pitchWeight * contribution;
    totalWeight += pitchWeight;
    nValidPoints++;
  }

  if (nValidPoints < 3) {
    return -9999.0;
  }

  if (totalWeight > 0) {
    logLikelihood = logLikelihood / totalWeight * nValidPoints;
  }

  return logLikelihood;
}

//***************************************************************
LikelihoodFitResult FitMomentumLikelihood(AnaParticlePD* particle,
                                          TProfile* dEdxTemplate,
                                          TGraph* rangeEnergyGraph,
                                          double particleMass) {
//***************************************************************

  LikelihoodFitResult result;
  result.momentum = -999;
  result.logLikelihood = -9999;
  result.uncertainty = -999;
  result.nHitsUsed = 0;
  result.valid = false;

  if (!particle || !dEdxTemplate || !rangeEnergyGraph) {
    return result;
  }

  const int plane = 2;
  if (particle->Hits[plane].empty() || particle->Hits[plane].size() < 3) {
    return result;
  }

  std::vector<double> measuredDeDx;
  std::vector<double> distanceFromVertex;
  std::vector<double> pitch;

  double cumulativeDistance = 0.0;
  TVector3 previousPos = particle->Hits[plane][0].Position;

  for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
    AnaHitPD& hit = particle->Hits[plane][i];

    if (hit.dEdx <= 0 || hit.dEdx > 1000 || hit.dEdx == -999) {
      continue;
    }

    if (i > 0) {
      if (hit.Pitch > 0 && hit.Pitch != -999) {
        cumulativeDistance += hit.Pitch;
      } else {
        TVector3 currentPos = hit.Position;
        cumulativeDistance += (currentPos - previousPos).Mag();
      }
      previousPos = hit.Position;
    }

    measuredDeDx.push_back(hit.dEdx);
    distanceFromVertex.push_back(cumulativeDistance);
    pitch.push_back(hit.Pitch);
  }

  if (measuredDeDx.size() < 3) {
    return result;
  }

  result.nHitsUsed = measuredDeDx.size();

  const double momentumMin = 0.05;
  const double momentumMax = 2.0;
  const double momentumStep = 0.01;

  double maxLogLikelihood = -9999.0;
  double bestMomentum = -999.0;

  for (double p = momentumMin; p <= momentumMax; p += momentumStep) {
    const double logL = CalculateLogLikelihood(measuredDeDx, distanceFromVertex, pitch, p,
                                               dEdxTemplate, rangeEnergyGraph, particleMass);
    if (logL > maxLogLikelihood) {
      maxLogLikelihood = logL;
      bestMomentum = p;
    }
  }

  result.momentum = bestMomentum;
  result.logLikelihood = maxLogLikelihood;

  double uncertaintyEstimate = -999.0;
  if (bestMomentum > 0) {
    const double targetLogL = maxLogLikelihood - 0.5;
    for (double dp = momentumStep; dp < 0.5; dp += momentumStep) {
      const double p_test = bestMomentum + dp;
      if (p_test > momentumMax) break;

      const double logL = CalculateLogLikelihood(measuredDeDx, distanceFromVertex, pitch, p_test,
                                                 dEdxTemplate, rangeEnergyGraph, particleMass);
      if (logL < targetLogL) {
        uncertaintyEstimate = dp;
        break;
      }
    }
  }
  result.uncertainty = uncertaintyEstimate;

  result.valid = (bestMomentum > 0 &&
                  bestMomentum < momentumMax &&
                  maxLogLikelihood > -9000 &&
                  result.nHitsUsed >= 3);

  return result;
}

//***************************************************************
Float_t EstimateMomentumByLikelihood(AnaParticlePD* particle, int pdg) {
//***************************************************************

  if (!particle || pdg != 211) {
    return -999.0;
  }

  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;
  pdMomShared::LoadPionTemplates(pionDEdxTemplate, pionRangeEnergyGraph,
                                 "pionRangeEnergyGraph_likelihood");

  if (!pionDEdxTemplate || !pionRangeEnergyGraph) {
    std::cout << "WARNING: Failed to load pion templates for likelihood momentum estimation" << std::endl;
    return -999.0;
  }

  LikelihoodFitResult result = FitMomentumLikelihood(particle, pionDEdxTemplate,
                                                     pionRangeEnergyGraph, kPionMassMeV);

  if (result.valid) {
    return static_cast<Float_t>(result.momentum);
  }

  return -999.0;
}

} // namespace pdMomReconstruction

namespace pdMomShared {

//***************************************************************
bool LoadPionTemplates(TProfile*& pionDEdxTemplate, TGraph*& pionRangeEnergyGraph,
                       const char* graphCloneName) {
//***************************************************************
  if (!pionDEdxTemplate) {
    TFile* templateFile = TFile::Open((std::string(getenv("PDUTILSROOT")) +
                                      "/data/dEdxrestemplates.root").c_str(), "READ");
    if (templateFile && !templateFile->IsZombie()) {
      pionDEdxTemplate = (TProfile*)templateFile->Get("dedx_range_pi");
      if (pionDEdxTemplate) {
        pionDEdxTemplate->SetDirectory(0);
      }
      templateFile->Close();
      delete templateFile;
    }
  }

  if (!pionRangeEnergyGraph) {
    TFile* rangeFile = TFile::Open((std::string(getenv("PDUTILSROOT")) +
                                   "/data/ke_vs_range.root").c_str(), "READ");
    if (rangeFile && !rangeFile->IsZombie()) {
      TGraph* tempGraph = (TGraph*)rangeFile->Get("pion");
      if (tempGraph) {
        pionRangeEnergyGraph = (TGraph*)tempGraph->Clone(graphCloneName);
      }
      rangeFile->Close();
      delete rangeFile;
    }
  }

  return (pionDEdxTemplate && pionRangeEnergyGraph);
}

//***************************************************************
double MomentumToKineticEnergyMeV(double momentumGeV, double massMeV) {
//***************************************************************
  const double momentumMeV = momentumGeV * 1000.0;
  const double totalEnergy = std::sqrt(momentumMeV * momentumMeV + massMeV * massMeV);
  return totalEnergy - massMeV;
}

//***************************************************************
double KineticEnergyMeVToRangeCm(double kineticEnergyMeV, TGraph* rangeEnergyGraph) {
//***************************************************************
  if (!rangeEnergyGraph || kineticEnergyMeV <= 0.) {
    return -999.0;
  }

  int nPoints = rangeEnergyGraph->GetN();
  double* xVals = rangeEnergyGraph->GetX();
  double* yVals = rangeEnergyGraph->GetY();

  for (int i = 0; i < nPoints - 1; i++) {
    if (yVals[i] <= kineticEnergyMeV && kineticEnergyMeV <= yVals[i + 1]) {
      const double slope = (xVals[i + 1] - xVals[i]) / (yVals[i + 1] - yVals[i]);
      return xVals[i] + slope * (kineticEnergyMeV - yVals[i]);
    }
  }

  if (kineticEnergyMeV > yVals[nPoints - 1]) {
    const double slope = (xVals[nPoints - 1] - xVals[nPoints - 2]) / (yVals[nPoints - 1] - yVals[nPoints - 2]);
    return xVals[nPoints - 1] + slope * (kineticEnergyMeV - yVals[nPoints - 1]);
  }

  return -999.0;
}

//***************************************************************
double RangeCmToKineticEnergyMeV(double rangeCm, TGraph* rangeEnergyGraph) {
//***************************************************************
  if (!rangeEnergyGraph || rangeCm < 0.) {
    return -999.0;
  }

  const double kineticEnergy = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(rangeEnergyGraph, rangeCm);
  if (kineticEnergy < 0. || kineticEnergy > 1e6 || !std::isfinite(kineticEnergy)) {
    return -999.0;
  }
  return kineticEnergy;
}

//***************************************************************
double RangeCmToMomentumGeV(double rangeCm, int /*pdg*/, TGraph* rangeEnergyGraph, double massMeV) {
//***************************************************************
  const double kineticEnergy = RangeCmToKineticEnergyMeV(rangeCm, rangeEnergyGraph);
  if (kineticEnergy < 0.) {
    return -999.0;
  }

  const double momentumMeV = std::sqrt(kineticEnergy * kineticEnergy + 2.0 * massMeV * kineticEnergy);
  return momentumMeV / 1000.0;
}

} // namespace pdMomShared

