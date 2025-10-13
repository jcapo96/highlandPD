#ifndef pdMomLikelihood_h
#define pdMomLikelihood_h

// Likelihood-based momentum estimation for non-stopping, interacting pions
// Specifically designed for K0 decay pions that don't exhibit Bragg peaks

#include "pdDataClasses.hxx"
#include <TProfile.h>
#include <TGraph.h>

namespace pdMomLikelihood {

  // Result structure for likelihood-based momentum estimation
  struct LikelihoodFitResult {
    double momentum;           // Best-fit momentum (GeV/c)
    double logLikelihood;      // Log-likelihood at best fit
    double uncertainty;        // Estimated uncertainty (GeV/c)
    int nHitsUsed;            // Number of hits used in fit
    bool valid;               // Whether fit was successful
  };

  // Main function: estimate momentum using likelihood method
  Float_t EstimateMomentumWithLikelihood(AnaParticlePD* particle, int pdg = 211);

  // Perform likelihood scan over momentum hypotheses
  LikelihoodFitResult FitMomentumLikelihood(AnaParticlePD* particle,
                                            TProfile* dEdxTemplate,
                                            TGraph* rangeEnergyGraph,
                                            double particleMass);

  // Calculate log-likelihood for a given momentum hypothesis
  double CalculateLogLikelihood(const std::vector<double>& measuredDeDx,
                                const std::vector<double>& distanceFromVertex,
                                const std::vector<double>& pitch,
                                double momentumHypothesis,
                                TProfile* dEdxVsRangeTemplate,
                                TGraph* rangeEnergyGraph,
                                double particleMass);

  // Get expected dE/dx at a given distance for a particle with given momentum
  double GetExpectedDeDx(double distanceTraveled,
                         double momentum,
                         TProfile* dEdxVsRangeTemplate,
                         TGraph* rangeEnergyGraph,
                         double particleMass);

  // Convert momentum to initial kinetic energy (MeV)
  double MomentumToKineticEnergy(double momentumGeV, double massMeV);

  // Calculate residual range from kinetic energy using range-energy table
  double KineticEnergyToRange(double kineticEnergyMeV, TGraph* rangeEnergyGraph);

  // Calculate kinetic energy from range using range-energy table
  double RangeToKineticEnergy(double rangeCm, TGraph* rangeEnergyGraph);

} // namespace pdMomLikelihood

#endif

