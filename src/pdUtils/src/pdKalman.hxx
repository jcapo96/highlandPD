#ifndef pdKalman_h
#define pdKalman_h

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <vector>
#include "pdDataClasses.hxx"

namespace pdKalman {

  // Track state at a point: position (3D) + direction (3D unit vector)
  struct TrackState {
    TVector3 position;
    TVector3 direction;
    TMatrixD covariance;  // 6x6 covariance matrix
  };

  // Vertex state: 3D position
  struct VertexState {
    TVector3 position;
    TMatrixD covariance;  // 3x3 covariance matrix
    double chi2;
    int ndf;
  };

  // Estimate track state covariance from fit residuals
  TMatrixD EstimateTrackCovariance(AnaParticlePD* particle, const TVector3& position,
                                   const TVector3& direction, double trackFitLength);

  // Kalman vertex fit for two tracks
  VertexState FitVertex(const TrackState& track1, const TrackState& track2,
                       const TVector3& seedPosition);

  // Convert AnaParticlePD to TrackState
  TrackState ParticleToTrackState(AnaParticlePD* particle, double trackFitLength,
                                  bool useStartPosition);
}

#endif

