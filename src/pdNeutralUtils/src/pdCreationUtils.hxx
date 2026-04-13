#ifndef pdCreationUtils_h
#define pdCreationUtils_h

#include "pdDataClasses.hxx"
#include <vector>
#include <TVector3.h>

namespace pdCreationUtils {

	TVector3 ProjectPointOntoLine(const TVector3& point,
																const TVector3& linePoint,
																const TVector3& lineDirection);

	bool ProjectBeamTailOntoStartDirection(const AnaParticlePD* beamParticle,
											 double fitDistanceFromEndCm,
																				 TVector3& correctedEndPosition,
																				 std::vector<TVector3>* rawTailPoints = nullptr,
																				 std::vector<TVector3>* projectedTailPoints = nullptr);
}

#endif