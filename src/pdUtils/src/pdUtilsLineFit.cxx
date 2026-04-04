#include "pdUtilsLineFit.hxx"
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TVectorD.h>
#include <algorithm>
#include <cmath>

//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition,
                                  double trackFitDistanceFromStart){

  fitParams.clear();
  fitParams.resize(6, -999.0);

  if (!part) {
    return;
  }

  TVector3 referencePos = pdAnaUtils::DefinePosition(part, useStartPosition);
  if (referencePos.X() < -900) {
    return;
  }

  TVector3 travelDir(-999, -999, -999);
  if (useStartPosition) {
    travelDir.SetXYZ(part->DirectionStart[0], part->DirectionStart[1], part->DirectionStart[2]);
  } else {
    travelDir.SetXYZ(-part->DirectionEnd[0], -part->DirectionEnd[1], -part->DirectionEnd[2]);
  }
  const bool hasValidTravelDir = (travelDir.Mag2() > 1e-10);
  if (hasValidTravelDir) travelDir = travelDir.Unit();

  std::vector<std::pair<TVector3, double>> trajectoryPointsWithDistance;
  if (part->TrjPoints.size() >= 2) {
    trajectoryPointsWithDistance.reserve(part->TrjPoints.size());
    double cumulative = 0.0;
    TVector3 prev;
    bool hasPrev = false;

    if (useStartPosition) {
      for (size_t i = 0; i < part->TrjPoints.size(); ++i) {
        const TVector3 pos = part->TrjPoints[i].Position;
        if (pos.Z() == -999) continue;
        if (hasPrev) cumulative += (pos - prev).Mag();
        trajectoryPointsWithDistance.push_back(std::make_pair(pos, cumulative));
        prev = pos;
        hasPrev = true;
      }
    } else {
      for (int i = static_cast<int>(part->TrjPoints.size()) - 1; i >= 0; --i) {
        const TVector3 pos = part->TrjPoints[i].Position;
        if (pos.Z() == -999) continue;
        if (hasPrev) cumulative += (pos - prev).Mag();
        trajectoryPointsWithDistance.push_back(std::make_pair(pos, cumulative));
        prev = pos;
        hasPrev = true;
      }
    }
  }

  auto computePathDistance = [&](const TVector3& position) -> double {
    if (!trajectoryPointsWithDistance.empty()) {
      double bestDist2 = 1e30;
      double bestArc = -1.0;
      for (const auto& tp : trajectoryPointsWithDistance) {
        const double d2 = (position - tp.first).Mag2();
        if (d2 < bestDist2) {
          bestDist2 = d2;
          bestArc = tp.second;
        }
      }
      return bestArc;
    }

    TVector3 delta = position - referencePos;
    return hasValidTravelDir ? delta.Dot(travelDir) : delta.Mag();
  };

  std::vector<std::pair<TVector3, double>> hitPositionsWithDistance;

  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part->Hits[plane].size(); i++) {
      AnaHitPD& hit = part->Hits[plane][i];
      TVector3 position;

      if (hit.Position.Z() != -999) {
        position = hit.Position;
      } else {
        continue;
      }

      const double pathDistance = computePathDistance(position);
      hitPositionsWithDistance.push_back(std::make_pair(position, pathDistance));
    }
  }

  if (hitPositionsWithDistance.size() < 2) {
    return;
  }

  const double fitWindowStart = std::max(0.0, trackFitDistanceFromStart);
  const double fitWindowEnd = fitWindowStart + std::max(0.0, trackLength);

  TVector3 anchorPoint = referencePos;
  bool foundAnchor = false;
  double bestAnchorDelta = 1e30;
  for (const auto& hitPair : hitPositionsWithDistance) {
    if (hitPair.second < 0.0) continue;
    const double delta = std::abs(hitPair.second - fitWindowStart);
    if (!foundAnchor || delta < bestAnchorDelta) {
      bestAnchorDelta = delta;
      anchorPoint = hitPair.first;
      foundAnchor = true;
    }
  }
  if (!foundAnchor && hasValidTravelDir) {
    anchorPoint = referencePos + fitWindowStart * travelDir;
  }

  std::vector<TVector3> nearbyHits;
  for (const auto& hitPair : hitPositionsWithDistance) {
    if (hasValidTravelDir) {
      if (hitPair.second >= fitWindowStart && hitPair.second <= fitWindowEnd) {
        nearbyHits.push_back(hitPair.first);
      }
    } else if (hitPair.second >= fitWindowStart && hitPair.second <= fitWindowEnd) {
      nearbyHits.push_back(hitPair.first);
    }
  }

  if (nearbyHits.size() >= 2) {
    TVector3 centroid(0, 0, 0);
    for (const auto& pos : nearbyHits) {
      centroid += pos;
    }
    centroid *= (1.0 / nearbyHits.size());

    double covXX = 0, covXY = 0, covXZ = 0;
    double covYY = 0, covYZ = 0, covZZ = 0;
    for (const auto& pos : nearbyHits) {
      const double dx = pos.X() - centroid.X();
      const double dy = pos.Y() - centroid.Y();
      const double dz = pos.Z() - centroid.Z();
      covXX += dx * dx;
      covXY += dx * dy;
      covXZ += dx * dz;
      covYY += dy * dy;
      covYZ += dy * dz;
      covZZ += dz * dz;
    }

    const int nPoints = nearbyHits.size();
    covXX /= nPoints;
    covXY /= nPoints;
    covXZ /= nPoints;
    covYY /= nPoints;
    covYZ /= nPoints;
    covZZ /= nPoints;

    TMatrixD covMatrix(3, 3);
    covMatrix(0, 0) = covXX; covMatrix(0, 1) = covXY; covMatrix(0, 2) = covXZ;
    covMatrix(1, 0) = covXY; covMatrix(1, 1) = covYY; covMatrix(1, 2) = covYZ;
    covMatrix(2, 0) = covXZ; covMatrix(2, 1) = covYZ; covMatrix(2, 2) = covZZ;

    TMatrixDEigen eigen(covMatrix);
    TVectorD eigenValues = eigen.GetEigenValuesRe();
    TMatrixD eigenVectors = eigen.GetEigenVectors();

    int maxEigenIndex = 0;
    double maxEigenValue = eigenValues[0];
    for (int i = 1; i < 3; i++) {
      if (eigenValues[i] > maxEigenValue) {
        maxEigenValue = eigenValues[i];
        maxEigenIndex = i;
      }
    }

    TVector3 direction(eigenVectors(0, maxEigenIndex),
                       eigenVectors(1, maxEigenIndex),
                       eigenVectors(2, maxEigenIndex));
    direction = direction.Unit();

    fitParams[0] = anchorPoint.X();
    fitParams[1] = anchorPoint.Y();
    fitParams[2] = anchorPoint.Z();
    fitParams[3] = direction.X();
    fitParams[4] = direction.Y();
    fitParams[5] = direction.Z();
    if (hasValidTravelDir) {
      TVector3 fitDir(fitParams[3], fitParams[4], fitParams[5]);
      if (fitDir.Mag2() > 1e-10 && fitDir.Dot(travelDir) < 0.0) {
        fitParams[3] *= -1.0;
        fitParams[4] *= -1.0;
        fitParams[5] *= -1.0;
      }
    }
  }

}

//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle, bool useStartPosition) {
//***************************************************************

  if (!particle) {
    return TVector3(-999, -999, -999);
  }

  if (useStartPosition) {
    return TVector3(particle->PositionStart[0],
                    particle->PositionStart[1],
                    particle->PositionStart[2]);
  } else {
    return TVector3(particle->PositionEnd[0],
                    particle->PositionEnd[1],
                    particle->PositionEnd[2]);
  }

}

//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle) {
//***************************************************************
  return DefinePosition(particle, true);
}

//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength) {
//***************************************************************
  ExtrapolateTrack(part, fitParams, trackLength, true, 0.0);
}

//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition){

  fitParams.clear();
  fitParams.resize(6, -999.0);

  if (!part) {
    return;
  }

  double x0, y0, z0, dx, dy, dz;

  if (useStartPosition) {
    x0 = part->Position[0];
    y0 = part->Position[1];
    z0 = part->Position[2];
    dx = part->Direction[0];
    dy = part->Direction[1];
    dz = part->Direction[2];
  } else {
    x0 = part->PositionEnd[0];
    y0 = part->PositionEnd[1];
    z0 = part->PositionEnd[2];
    dx = part->DirectionEnd[0];
    dy = part->DirectionEnd[1];
    dz = part->DirectionEnd[2];
  }

  if (x0 == -999.0 || y0 == -999.0 || z0 == -999.0 ||
      dx == -999.0 || dy == -999.0 || dz == -999.0) {
    return;
  }

  double norm = sqrt(dx*dx + dy*dy + dz*dz);
  if (norm > 0) {
    dx /= norm;
    dy /= norm;
    dz /= norm;
  } else {
    return;
  }

  fitParams[0] = x0;
  fitParams[1] = y0;
  fitParams[2] = z0;
  fitParams[3] = dx;
  fitParams[4] = dy;
  fitParams[5] = dz;

}

//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength) {
//***************************************************************
  ExtrapolateTrack(part, fitParams, trackLength, true);
}

//***************************************************************
double pdAnaUtils::FindClosestPointsBetweenLines(const std::vector<double>& line1Params,
                                                const std::vector<double>& line2Params,
                                                TVector3& closestPoint1,
                                                TVector3& closestPoint2) {

  if (line1Params.size() != 6 || line2Params.size() != 6 ||
      line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    closestPoint1.SetXYZ(-999, -999, -999);
    closestPoint2.SetXYZ(-999, -999, -999);
    return -999.0;
  }

  TVector3 point1(line1Params[0], line1Params[1], line1Params[2]);
  TVector3 dir1(line1Params[3], line1Params[4], line1Params[5]);

  TVector3 point2(line2Params[0], line2Params[1], line2Params[2]);
  TVector3 dir2(line2Params[3], line2Params[4], line2Params[5]);

  dir1 = dir1.Unit();
  dir2 = dir2.Unit();

  TVector3 w0 = point1 - point2;

  double a = dir1.Dot(dir1);
  double b = dir1.Dot(dir2);
  double c = dir2.Dot(dir2);
  double d = dir1.Dot(w0);
  double e = dir2.Dot(w0);

  double denom = a * c - b * b;

  if (fabs(denom) < 1e-10) {
    closestPoint1 = point1;
    double t2 = e / c;
    closestPoint2 = point2 + t2 * dir2;
  } else {
    double t1 = (b * e - c * d) / denom;
    double t2 = (a * e - b * d) / denom;

    closestPoint1 = point1 + t1 * dir1;
    closestPoint2 = point2 + t2 * dir2;
  }

  double minDistance = (closestPoint1 - closestPoint2).Mag();

  return minDistance;
}

//********************************************************************
double pdAnaUtils::CalculateImpactParameter(const std::vector<double>& lineParams, const TVector3& point){
//********************************************************************

  if (lineParams.size() != 6 || lineParams[0] == -999.0) {
    return -999.0;
  }

  TVector3 linePoint(lineParams[0], lineParams[1], lineParams[2]);
  TVector3 lineDirection(lineParams[3], lineParams[4], lineParams[5]);
  lineDirection = lineDirection.Unit();

  TVector3 pointToLine = point - linePoint;
  TVector3 projection = (pointToLine.Dot(lineDirection)) * lineDirection;
  TVector3 perpendicular = pointToLine - projection;

  return perpendicular.Mag();
}
