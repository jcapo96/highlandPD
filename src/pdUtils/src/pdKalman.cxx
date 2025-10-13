#include "pdKalman.hxx"
#include "pdAnalysisUtils.hxx"
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMath.h>
#include <iostream>

namespace pdKalman {

//***************************************************************
TMatrixD EstimateTrackCovariance(AnaParticlePD* particle, const TVector3& position,
                                 const TVector3& direction, double trackFitLength) {
//***************************************************************

  // Create and properly initialize 6x6 covariance matrix (position + direction)
  TMatrixD covariance(6, 6);

  // Initialize to zero first
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      covariance(i, j) = 0.0;
    }
  }

  if (!particle) {
    // Return large diagonal covariance if no particle
    for (int i = 0; i < 3; i++) {
      covariance(i, i) = 100.0;  // 10 cm position uncertainty
      covariance(i+3, i+3) = 0.01;  // 0.1 rad direction uncertainty
    }
    return covariance;
  }

  // Collect hits within trackFitLength of reference position
  std::vector<TVector3> nearbyHits;

  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < particle->Hits[plane].size(); i++) {
      AnaHitPD& hit = particle->Hits[plane][i];
      if (hit.Position.Z() != -999) {
        TVector3 hitPos = hit.Position;
        double distance = (hitPos - position).Mag();
        if (distance <= trackFitLength) {
          nearbyHits.push_back(hitPos);
        }
      }
    }
  }

  // Need at least 3 hits to estimate uncertainties
  if (nearbyHits.size() < 3) {
    // Use default uncertainties
    for (int i = 0; i < 3; i++) {
      covariance(i, i) = 4.0;  // 2 cm position uncertainty
      covariance(i+3, i+3) = 0.0004;  // 0.02 rad direction uncertainty
    }
    return covariance;
  }

  // Calculate position uncertainties from hit scatter
  TVector3 meanPos(0, 0, 0);
  for (const auto& hit : nearbyHits) {
    meanPos += hit;
  }
  meanPos *= (1.0 / nearbyHits.size());

  // Calculate RMS in each direction
  double rmsX = 0, rmsY = 0, rmsZ = 0;
  for (const auto& hit : nearbyHits) {
    rmsX += pow(hit.X() - meanPos.X(), 2);
    rmsY += pow(hit.Y() - meanPos.Y(), 2);
    rmsZ += pow(hit.Z() - meanPos.Z(), 2);
  }
  rmsX = sqrt(rmsX / nearbyHits.size());
  rmsY = sqrt(rmsY / nearbyHits.size());
  rmsZ = sqrt(rmsZ / nearbyHits.size());

  // Add systematic uncertainty
  double syst = 0.5; // cm
  rmsX = sqrt(rmsX*rmsX + syst*syst);
  rmsY = sqrt(rmsY*rmsY + syst*syst);
  rmsZ = sqrt(rmsZ*rmsZ + syst*syst);

  // Ensure minimum uncertainty
  if (rmsX < 0.5) rmsX = 0.5;
  if (rmsY < 0.5) rmsY = 0.5;
  if (rmsZ < 0.5) rmsZ = 0.5;

  // Set position covariance (variance = rms^2)
  covariance(0, 0) = rmsX * rmsX;
  covariance(1, 1) = rmsY * rmsY;
  covariance(2, 2) = rmsZ * rmsZ;

  // Estimate direction uncertainties from angular spread
  // Calculate perpendicular distances from hits to line
  std::vector<double> perpDistances;
  for (const auto& hit : nearbyHits) {
    TVector3 toHit = hit - position;
    TVector3 projection = direction * (toHit.Dot(direction));
    TVector3 perpendicular = toHit - projection;
    perpDistances.push_back(perpendicular.Mag());
  }

  // RMS perpendicular distance
  double rmsPerpDist = 0;
  for (double d : perpDistances) {
    rmsPerpDist += d * d;
  }
  rmsPerpDist = sqrt(rmsPerpDist / perpDistances.size());

  // Angular uncertainty ~ perpDist / trackLength
  double angularUncertainty = rmsPerpDist / trackFitLength;
  if (angularUncertainty < 0.01) angularUncertainty = 0.01;  // Min 0.01 rad
  if (angularUncertainty > 0.1) angularUncertainty = 0.1;    // Max 0.1 rad

  // Set direction covariance (isotropic in angle space)
  double dirVariance = angularUncertainty * angularUncertainty;
  covariance(3, 3) = dirVariance;
  covariance(4, 4) = dirVariance;
  covariance(5, 5) = dirVariance;

  return covariance;
}

//***************************************************************
VertexState FitVertex(const TrackState& track1, const TrackState& track2,
                     const TVector3& seedPosition) {
//***************************************************************

  // Following Fruhwirth's Kalman filter vertex fit formalism
  // Nucl. Inst. Meth. A262 (1987) 444

  VertexState result;

  // Initialize vertex state with seed position and large covariance
  TVectorD x_vtx(3);  // Vertex position
  x_vtx(0) = seedPosition.X();
  x_vtx(1) = seedPosition.Y();
  x_vtx(2) = seedPosition.Z();

  TMatrixD C_vtx(3, 3);  // Vertex covariance
  C_vtx.Zero();
  C_vtx(0, 0) = 100.0;  // 10 cm initial uncertainty
  C_vtx(1, 1) = 100.0;
  C_vtx(2, 2) = 100.0;

  double totalChi2 = 0.0;
  int ndf = 0;

  // Process track 1
  {
    // Measurement: track position (extracted from track state)
    TVectorD z1(3);
    z1(0) = track1.position.X();
    z1(1) = track1.position.Y();
    z1(2) = track1.position.Z();

    // Measurement covariance: position part of track covariance (3x3)
    TMatrixD V1(3, 3);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        V1(i, j) = track1.covariance(i, j);
      }
    }

    // Measurement model: H = I (identity, vertex position directly measured)
    TMatrixD H(3, 3);
    H.UnitMatrix();

    // Innovation covariance: S = H * C * H^T + V
    TMatrixD S = C_vtx + V1;

    // Kalman gain: K = C * H^T * S^-1
    TMatrixD S_inv = S;
    S_inv.Invert();
    TMatrixD K = C_vtx * S_inv;

    // Residual: r = z - x_pred
    TVectorD residual = z1 - x_vtx;

    // Update state: x = x + K * r
    x_vtx += K * residual;

    // Update covariance: C = (I - K * H) * C = (I - K) * C
    TMatrixD I(3, 3);
    I.UnitMatrix();
    TMatrixD IminusK = I - K;
    C_vtx = IminusK * C_vtx;

    // Chi-squared: chi2 = r^T * S^-1 * r
    TVectorD temp = S_inv * residual;
    double chi2_contrib = residual * temp;
    totalChi2 += chi2_contrib;
    ndf += 3;  // 3 measurements from this track
  }

  // Process track 2
  {
    // Measurement: track position
    TVectorD z2(3);
    z2(0) = track2.position.X();
    z2(1) = track2.position.Y();
    z2(2) = track2.position.Z();

    // Measurement covariance
    TMatrixD V2(3, 3);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        V2(i, j) = track2.covariance(i, j);
      }
    }

    // Measurement model: H = I
    TMatrixD H(3, 3);
    H.UnitMatrix();

    // Innovation covariance
    TMatrixD S = C_vtx + V2;

    // Kalman gain
    TMatrixD S_inv = S;
    S_inv.Invert();
    TMatrixD K = C_vtx * S_inv;

    // Residual
    TVectorD residual = z2 - x_vtx;

    // Update state
    x_vtx += K * residual;

    // Update covariance
    TMatrixD I(3, 3);
    I.UnitMatrix();
    TMatrixD IminusK = I - K;
    C_vtx = IminusK * C_vtx;

    // Chi-squared
    TVectorD temp = S_inv * residual;
    double chi2_contrib = residual * temp;
    totalChi2 += chi2_contrib;
    ndf += 3;  // 3 measurements from this track
  }

  // Degrees of freedom: 2 tracks * 3 measurements - 3 vertex parameters
  ndf -= 3;

  // Store results
  result.position.SetXYZ(x_vtx(0), x_vtx(1), x_vtx(2));
  result.covariance.ResizeTo(3, 3);
  result.covariance = C_vtx;
  result.chi2 = totalChi2;
  result.ndf = ndf;

  return result;
}

//***************************************************************
TrackState ParticleToTrackState(AnaParticlePD* particle, double trackFitLength,
                                bool useStartPosition) {
//***************************************************************

  TrackState state;

  // Get position and direction from particle
  if (!particle) {
    state.position.SetXYZ(-999, -999, -999);
    state.direction.SetXYZ(0, 0, 1);
  } else if (useStartPosition) {
    state.position.SetXYZ(particle->PositionStart[0],
                         particle->PositionStart[1],
                         particle->PositionStart[2]);
    state.direction.SetXYZ(particle->DirectionStart[0],
                          particle->DirectionStart[1],
                          particle->DirectionStart[2]);
  } else {
    state.position.SetXYZ(particle->PositionEnd[0],
                         particle->PositionEnd[1],
                         particle->PositionEnd[2]);
    state.direction.SetXYZ(particle->DirectionEnd[0],
                          particle->DirectionEnd[1],
                          particle->DirectionEnd[2]);
  }

  // Normalize direction
  if (state.direction.Mag() > 0) {
    state.direction = state.direction.Unit();
  } else {
    // Default to z direction if invalid
    state.direction.SetXYZ(0, 0, 1);
  }

  // Estimate covariance matrix and resize state.covariance before assignment
  TMatrixD estimatedCov = EstimateTrackCovariance(particle, state.position,
                                                  state.direction, trackFitLength);
  state.covariance.ResizeTo(6, 6);
  state.covariance = estimatedCov;

  return state;
}

} // namespace pdKalman

