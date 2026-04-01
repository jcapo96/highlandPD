#include "pdMomEstimation.hxx"
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TMath.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// Track-length extension fitting algorithm for momentum estimation of non-stopping particles
// Based on: "The track-length extension fitting algorithm for energy measurement
// of interacting particles in liquid argon TPCs and its performance with ProtoDUNE-SP data"
// DUNE Collaboration, arXiv:2409.18288
// https://arxiv.org/pdf/2409.18288

namespace pdMomEstimation {

namespace {
bool LoadPionTemplates(TProfile*& pionDEdxTemplate, TGraph*& pionRangeEnergyGraph) {
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
        pionRangeEnergyGraph = (TGraph*)tempGraph->Clone("pionRangeEnergyGraph_clone");
      }
      rangeFile->Close();
      delete rangeFile;
    }
  }

  return (pionDEdxTemplate && pionRangeEnergyGraph);
}
} // namespace

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

  // Calculate chi2 comparing measured dE/dx to template with extension offset
  for (size_t i = 0; i < measuredDeDx.size(); i++) {
    // Effective residual range = measured RR + extension
    double effectiveRR = measuredRR[i] + extension;

    // Get expected dE/dx from template at effective RR
    double expectedDeDx = dEdxTemplate->Interpolate(effectiveRR);

    // Skip if template value is invalid or outside range
    if (expectedDeDx <= 0 || expectedDeDx > 1000) {
      continue;
    }

    // Measurement uncertainty (same model as Chi2PID)
    double sigma = 0.04231 + 0.0001783 * measuredDeDx[i] * measuredDeDx[i];
    sigma *= measuredDeDx[i];

    // Add minimum uncertainty to avoid division by zero
    if (sigma < 0.1) sigma = 0.1;

    // Chi-squared contribution
    double residual = measuredDeDx[i] - expectedDeDx;
    chi2 += (residual * residual) / (sigma * sigma);
    nPoints++;
  }

  // Penalize if too few points
  if (nPoints < 3) {
    return 9999.0;
  }

  return chi2;
}

//***************************************************************
double RangeToMomentum(double effectiveRange, int pdg, TGraph* rangeEnergyGraph, double mass) {
//***************************************************************

  if (!rangeEnergyGraph || effectiveRange <= 0) {
    return -999.0;
  }

  // Get kinetic energy from range using the graph
  // Note: graph should have range on x-axis, KE on y-axis
  double kineticEnergy = rangeEnergyGraph->Eval(effectiveRange);

  if (kineticEnergy <= 0 || kineticEnergy > 1e6) {
    return -999.0;
  }

  // Calculate momentum from kinetic energy and mass
  // KE is in MeV, mass is in MeV
  // p = sqrt(KE^2 + 2*m*KE)
  double momentum = sqrt(kineticEnergy * kineticEnergy + 2.0 * mass * kineticEnergy);

  // Convert to GeV/c
  momentum = momentum / 1000.0;

  return momentum;
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

  int plane = 2; // Collection plane

  // Check if particle has hits
  if (particle->Hits[plane].empty() || particle->Hits[plane].size() < 3) {
    return result;
  }

  // Extract measured dE/dx and calculate residual range
  std::vector<double> measuredDeDx;
  std::vector<double> distanceFromStart; // Distance from vertex along track

  // First pass: calculate distances from start and collect valid hits
  double cumulativeDistance = 0.0;
  TVector3 previousPos = particle->Hits[plane][0].Position;

  for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
    AnaHitPD& hit = particle->Hits[plane][i];

    // Skip invalid hits
    if (hit.dEdx <= 0 || hit.dEdx > 1000 || hit.dEdx == -999) {
      continue;
    }

    // Calculate distance from previous hit
    if (i > 0) {
      TVector3 currentPos = hit.Position;
      cumulativeDistance += (currentPos - previousPos).Mag();
      previousPos = currentPos;
    }

    measuredDeDx.push_back(hit.dEdx);
    distanceFromStart.push_back(cumulativeDistance);
  }

  // Second pass: calculate CORRECT residual range
  // RR = distance remaining to END of track
  // For stopping particles: RR=0 at end (Bragg peak), RR=max at start
  // For non-stopping particles: we extend beyond the measured track length
  std::vector<double> measuredRR;
  double totalMeasuredLength = distanceFromStart.empty() ? 0.0 : distanceFromStart.back();

  for (size_t i = 0; i < distanceFromStart.size(); i++) {
    // Residual range = total length - distance from start
    measuredRR.push_back(totalMeasuredLength - distanceFromStart[i]);
  }

  // Need at least 3 valid points
  if (measuredDeDx.size() < 3) {
    return result;
  }

  // Grid search for best extension
  double minChi2 = 9999.0;
  double bestExtension = 0.0;

  // Search range: 0 to 200 cm with 0.5 cm steps
  double extensionMin = 0.0;
  double extensionMax = 200.0;
  double extensionStep = 0.5;

  for (double ext = extensionMin; ext <= extensionMax; ext += extensionStep) {
    double chi2 = CalculateExtensionChi2(measuredDeDx, measuredRR, ext, dEdxTemplate);

    if (chi2 < minChi2) {
      minChi2 = chi2;
      bestExtension = ext;
    }
  }

  // Calculate effective range
  double measuredLength = measuredRR.back(); // Last cumulative distance
  double effectiveRange = measuredLength + bestExtension;

  // Calculate momentum from effective range
  double momentum = RangeToMomentum(effectiveRange, 211, rangeEnergyGraph, pionMass);

  // Store results
  result.extension = bestExtension;
  result.effectiveRange = effectiveRange;
  result.momentum = momentum;
  result.chi2 = minChi2;
  result.ndf = static_cast<int>(measuredDeDx.size());

  // Check if fit is reasonable
  double chi2ndf = (result.ndf > 0) ? result.chi2 / result.ndf : 9999.0;
  result.valid = (chi2ndf < 5.0 && momentum > 0 && momentum < 10.0); // Reasonable momentum range

  return result;
}

//***************************************************************
Float_t EstimateMomentumWithExtension(AnaParticlePD* particle, int pdg) {
//***************************************************************
  ExtensionFitResult result = EstimateMomentumWithExtensionDetailed(particle, pdg);
  if (result.valid) {
    return static_cast<Float_t>(result.momentum);
  }
  return -999.0;
}

//***************************************************************
Float_t EstimateMomentumFromPionRange(AnaParticlePD* particle) {
//***************************************************************
  if (!particle) return -999.0f;
  if (!(particle->Length > 0.0f) || particle->Length == -999.0f) return -999.0f;

  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;
  static constexpr double pionMass = 139.57; // MeV/c^2
  if (!LoadPionTemplates(pionDEdxTemplate, pionRangeEnergyGraph)) {
    return -999.0f;
  }

  const double momentum = RangeToMomentum(static_cast<double>(particle->Length), 211,
                                          pionRangeEnergyGraph, pionMass);
  if (!(momentum > 0.0) || !std::isfinite(momentum)) return -999.0f;
  return static_cast<Float_t>(momentum);
}

//***************************************************************
ExtensionFitResult EstimateMomentumWithExtensionDetailed(AnaParticlePD* particle, int pdg) {
//***************************************************************
  ExtensionFitResult invalidResult;
  invalidResult.extension = -999;
  invalidResult.effectiveRange = -999;
  invalidResult.momentum = -999;
  invalidResult.chi2 = 9999;
  invalidResult.ndf = 0;
  invalidResult.valid = false;

  if (!particle || pdg != 211) {
    // Only implemented for pions currently
    return invalidResult;
  }

  static TProfile* pionDEdxTemplate = nullptr;
  static TGraph* pionRangeEnergyGraph = nullptr;
  if (!LoadPionTemplates(pionDEdxTemplate, pionRangeEnergyGraph)) {
    std::cout << "WARNING: Failed to load pion templates for momentum estimation" << std::endl;
    return invalidResult;
  }

  const double pionMass = 139.57;
  return FitTrackLengthExtension(particle, pionDEdxTemplate, pionRangeEnergyGraph, pionMass);
}

} // namespace pdMomEstimation

