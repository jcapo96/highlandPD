#ifndef pdMomEstimation_h
#define pdMomEstimation_h

// Track-length extension fitting algorithm for momentum estimation of non-stopping particles
// Based on: "The track-length extension fitting algorithm for energy measurement
// of interacting particles in liquid argon TPCs and its performance with ProtoDUNE-SP data"
// DUNE Collaboration, arXiv:2409.18288
// https://arxiv.org/pdf/2409.18288

#include "pdDataClasses.hxx"
#include <TProfile.h>
#include <TGraph.h>

namespace pdMomEstimation {

  // Result structure for track-length extension fit
  struct ExtensionFitResult {
    double extension;        // Track extension length (cm)
    double effectiveRange;   // Total effective range (cm)
    double momentum;         // Estimated momentum (GeV/c)
    double chi2;            // Chi-squared of fit
    int ndf;                // Number of degrees of freedom
    bool valid;             // Whether fit was successful
  };

  // Estimate momentum using track-length extension method
  Float_t EstimateMomentumWithExtension(AnaParticlePD* particle, int pdg = 211);

  // Perform track-length extension fit
  ExtensionFitResult FitTrackLengthExtension(AnaParticlePD* particle,
                                             TProfile* dEdxTemplate,
                                             TGraph* rangeEnergyGraph,
                                             double pionMass);

  // Calculate chi-squared for a given extension
  double CalculateExtensionChi2(const std::vector<double>& measuredDeDx,
                               const std::vector<double>& measuredRR,
                               double extension,
                               TProfile* dEdxTemplate);

  // Convert effective range to momentum using range-energy tables
  double RangeToMomentum(double effectiveRange, int pdg, TGraph* rangeEnergyGraph, double mass);
}

#endif

