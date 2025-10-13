#include "pdMomLikelihood.hxx"
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// Likelihood-based momentum estimation for non-stopping, interacting pions
// Designed for K0 decay pions that don't stop and don't show Bragg peak

namespace pdMomLikelihood {

//***************************************************************
double MomentumToKineticEnergy(double momentumGeV, double massMeV) {
//***************************************************************
  // Convert momentum (GeV/c) to kinetic energy (MeV)
  // KE = sqrt(p²c² + m²c⁴) - mc²
  double momentumMeV = momentumGeV * 1000.0; // Convert to MeV/c
  double totalEnergy = sqrt(momentumMeV*momentumMeV + massMeV*massMeV);
  double kineticEnergy = totalEnergy - massMeV;
  return kineticEnergy;
}

//***************************************************************
double KineticEnergyToRange(double kineticEnergyMeV, TGraph* rangeEnergyGraph) {
//***************************************************************
  if (!rangeEnergyGraph || kineticEnergyMeV <= 0) {
    return -999.0;
  }

  // Invert the range-energy relationship
  // Graph has range (cm) on x-axis, KE (MeV) on y-axis
  // We need to find range given KE

  int nPoints = rangeEnergyGraph->GetN();
  double* xVals = rangeEnergyGraph->GetX();
  double* yVals = rangeEnergyGraph->GetY();

  // Find bracketing points
  for (int i = 0; i < nPoints - 1; i++) {
    if (yVals[i] <= kineticEnergyMeV && kineticEnergyMeV <= yVals[i+1]) {
      // Linear interpolation
      double slope = (xVals[i+1] - xVals[i]) / (yVals[i+1] - yVals[i]);
      double range = xVals[i] + slope * (kineticEnergyMeV - yVals[i]);
      return range;
    }
  }

  // Extrapolate if outside range
  if (kineticEnergyMeV > yVals[nPoints-1]) {
    // Extrapolate linearly from last two points
    double slope = (xVals[nPoints-1] - xVals[nPoints-2]) / (yVals[nPoints-1] - yVals[nPoints-2]);
    double range = xVals[nPoints-1] + slope * (kineticEnergyMeV - yVals[nPoints-1]);
    return range;
  }

  return -999.0;
}

//***************************************************************
double RangeToKineticEnergy(double rangeCm, TGraph* rangeEnergyGraph) {
//***************************************************************
  if (!rangeEnergyGraph || rangeCm <= 0) {
    return -999.0;
  }

  double kineticEnergy = rangeEnergyGraph->Eval(rangeCm);

  if (kineticEnergy <= 0 || kineticEnergy > 1e6) {
    return -999.0;
  }

  return kineticEnergy;
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

  // Calculate initial kinetic energy from momentum
  double initialKE = MomentumToKineticEnergy(momentumGeV, particleMass);
  if (initialKE <= 0) {
    return -999.0;
  }

  // For short distances, use simple approximation
  if (distanceTraveled < 1.0) {
    double totalRange = KineticEnergyToRange(initialKE, rangeEnergyGraph);
    if (totalRange <= 0) return -999.0;
    double residualRange = totalRange - distanceTraveled;
    if (residualRange < 0) return -999.0;
    return dEdxVsRangeTemplate->Interpolate(residualRange);
  }

  // For longer distances, integrate energy loss step-by-step
  // This accounts for dE/dx variation as particle slows down
  double energyLost = 0.0;
  double stepSize = 0.5; // cm
  double currentKE = initialKE;

  // Integrate from start to distanceTraveled
  for (double d = 0; d < distanceTraveled; d += stepSize) {
    // Calculate current kinetic energy
    currentKE = initialKE - energyLost;

    // Check if particle would have stopped
    if (currentKE <= 0) {
      return -999.0; // Particle stopped before reaching this distance
    }

    // Get current residual range
    double currentRange = KineticEnergyToRange(currentKE, rangeEnergyGraph);
    if (currentRange <= 0) {
      return -999.0;
    }

    // Get dE/dx at current position
    double currentDeDx = dEdxVsRangeTemplate->Interpolate(currentRange);
    if (currentDeDx <= 0 || currentDeDx > 1000) {
      return -999.0;
    }

    // Accumulate energy loss: dE = dE/dx * dx
    energyLost += currentDeDx * stepSize;
  }

  // Return expected dE/dx at the final position
  double finalKE = initialKE - energyLost;
  if (finalKE <= 0) {
    return -999.0;
  }

  double finalRange = KineticEnergyToRange(finalKE, rangeEnergyGraph);
  if (finalRange <= 0) {
    return -999.0;
  }

  double finalDeDx = dEdxVsRangeTemplate->Interpolate(finalRange);

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

  // Nominal pitch for collection plane (cm)
  const double nominalPitch = 0.4792;

  // Calculate likelihood for each measurement
  for (size_t i = 0; i < measuredDeDx.size(); i++) {
    double measured = measuredDeDx[i];
    double distance = distanceFromVertex[i];
    double hitPitch = pitch[i];

    // Get expected dE/dx for this momentum hypothesis
    double expected = GetExpectedDeDx(distance, momentumHypothesis,
                                      dEdxVsRangeTemplate, rangeEnergyGraph,
                                      particleMass);

    // Skip if prediction failed
    if (expected <= 0 || expected == -999.0) {
      continue;
    }

    // Measurement uncertainty (same model as Chi2PID)
    double sigma = 0.04231 + 0.0001783 * measured * measured;
    sigma *= measured;

    // Add minimum uncertainty
    if (sigma < 0.1) sigma = 0.1;

    // Calculate pitch-based weight
    // Hits with pitch close to nominal are more reliable
    double pitchWeight = 1.0;
    if (hitPitch > 0 && hitPitch != -999) {
      double pitchDeviation = fabs(hitPitch - nominalPitch) / nominalPitch;
      pitchWeight = 1.0 / (1.0 + pitchDeviation);
    }

    // Gaussian log-likelihood with pitch weighting
    double residual = measured - expected;
    double contribution = -0.5 * (residual * residual) / (sigma * sigma);
    logLikelihood += pitchWeight * contribution;
    totalWeight += pitchWeight;
    nValidPoints++;
  }

  // Penalize if too few valid points
  if (nValidPoints < 3) {
    return -9999.0;
  }

  // Normalize by total weight to keep likelihood values comparable
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

  int plane = 2; // Collection plane

  // Check if particle has hits
  if (particle->Hits[plane].empty() || particle->Hits[plane].size() < 3) {
    return result;
  }

  // Extract measured dE/dx, distance from vertex, and pitch
  std::vector<double> measuredDeDx;
  std::vector<double> distanceFromVertex;
  std::vector<double> pitch;

  double cumulativeDistance = 0.0;
  TVector3 previousPos = particle->Hits[plane][0].Position;

  for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
    AnaHitPD& hit = particle->Hits[plane][i];

    // Skip invalid hits
    if (hit.dEdx <= 0 || hit.dEdx > 1000 || hit.dEdx == -999) {
      continue;
    }

    // Calculate distance from vertex using pitch (true path length through LAr)
    if (i > 0) {
      // Use pitch as the effective path length for this hit segment
      // Pitch represents the actual distance traveled through the detector
      if (hit.Pitch > 0 && hit.Pitch != -999) {
        cumulativeDistance += hit.Pitch;
      } else {
        // Fallback to 3D Euclidean distance if pitch not available
        TVector3 currentPos = hit.Position;
        cumulativeDistance += (currentPos - previousPos).Mag();
      }
      previousPos = hit.Position;
    }

    measuredDeDx.push_back(hit.dEdx);
    distanceFromVertex.push_back(cumulativeDistance);
    pitch.push_back(hit.Pitch);
  }

  // Need at least 3 valid points
  if (measuredDeDx.size() < 3) {
    return result;
  }

  result.nHitsUsed = measuredDeDx.size();

  // Scan momentum hypotheses to find maximum likelihood
  double momentumMin = 0.05;  // GeV/c (50 MeV/c)
  double momentumMax = 2.0;   // GeV/c (2000 MeV/c)
  double momentumStep = 0.01; // GeV/c (10 MeV/c steps)

  double maxLogLikelihood = -9999.0;
  double bestMomentum = -999.0;

  for (double p = momentumMin; p <= momentumMax; p += momentumStep) {
    double logL = CalculateLogLikelihood(measuredDeDx, distanceFromVertex, pitch, p,
                                        dEdxTemplate, rangeEnergyGraph, particleMass);

    if (logL > maxLogLikelihood) {
      maxLogLikelihood = logL;
      bestMomentum = p;
    }
  }

  // Store results
  result.momentum = bestMomentum;
  result.logLikelihood = maxLogLikelihood;

  // Estimate uncertainty from likelihood curvature
  // Scan near maximum to find where log-likelihood drops by 0.5 (1-sigma)
  double uncertaintyEstimate = -999.0;
  if (bestMomentum > 0) {
    double targetLogL = maxLogLikelihood - 0.5;

    // Scan upward
    for (double dp = momentumStep; dp < 0.5; dp += momentumStep) {
      double p_test = bestMomentum + dp;
      if (p_test > momentumMax) break;

      double logL = CalculateLogLikelihood(measuredDeDx, distanceFromVertex, pitch, p_test,
                                          dEdxTemplate, rangeEnergyGraph, particleMass);
      if (logL < targetLogL) {
        uncertaintyEstimate = dp;
        break;
      }
    }
  }
  result.uncertainty = uncertaintyEstimate;

  // Validate result
  result.valid = (bestMomentum > 0 &&
                  bestMomentum < momentumMax &&
                  maxLogLikelihood > -9000 &&
                  result.nHitsUsed >= 3);

  return result;
}

//***************************************************************
Float_t EstimateMomentumWithLikelihood(AnaParticlePD* particle, int pdg) {
//***************************************************************

  if (!particle || pdg != 211) {
    // Only implemented for pions currently
    return -999.0;
  }

  // Load pion dE/dx template and range-energy graph
  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;

  if (!pionDEdxTemplate) {
    TFile* templateFile = TFile::Open((std::string(getenv("PDUTILSROOT")) +
                                      "/data/dEdxrestemplates.root").c_str(), "READ");
    if (templateFile && !templateFile->IsZombie()) {
      pionDEdxTemplate = (TProfile*)templateFile->Get("dedx_range_pi");
      if (pionDEdxTemplate) {
        pionDEdxTemplate->SetDirectory(0); // Detach from file
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
        // Clone the graph so it persists after file is closed
        pionRangeEnergyGraph = (TGraph*)tempGraph->Clone("pionRangeEnergyGraph_likelihood");
      }
      rangeFile->Close();
      delete rangeFile;
    }
  }

  // Check if templates loaded successfully
  if (!pionDEdxTemplate || !pionRangeEnergyGraph) {
    std::cout << "WARNING: Failed to load pion templates for likelihood momentum estimation" << std::endl;
    return -999.0;
  }

  // Pion mass in MeV
  const double pionMass = 139.57;

  // Perform likelihood-based momentum fit
  LikelihoodFitResult result = FitMomentumLikelihood(particle, pionDEdxTemplate,
                                                     pionRangeEnergyGraph, pionMass);

  // Return momentum if fit was successful
  if (result.valid) {
    return static_cast<Float_t>(result.momentum);
  }

  return -999.0;
}

} // namespace pdMomLikelihood

