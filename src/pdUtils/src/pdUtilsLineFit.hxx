#ifndef pdUtilsLineFit_h
#define pdUtilsLineFit_h

#include "pdDataClasses.hxx"
#include <vector>

namespace pdAnaUtils {

  TVector3 DefinePosition(AnaParticlePD* particle, bool useStartPosition);

  TVector3 DefinePosition(AnaParticlePD* particle);

  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition,
                        double trackFitDistanceFromStart = 0.0);

  void ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength);

  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition);

  void ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength);

  double FindClosestPointsBetweenLines(const std::vector<double>& line1Params,
                                       const std::vector<double>& line2Params,
                                       TVector3& closestPoint1,
                                       TVector3& closestPoint2);

  double CalculateImpactParameter(const std::vector<double>& lineParams, const TVector3& point);

}

#endif
