#ifndef pdNeutralHelpers_h
#define pdNeutralHelpers_h

#include "pdDataClasses.hxx"
#include "TVector3.h"

namespace pdNeutralHelpers {

  // Calculate midpoint of minimum distance between two 3D lines
  // Line 1: P1(s) = pos1 + s * dir1
  // Line 2: P2(t) = pos2 + t * dir2
  // Returns the midpoint between the two closest points on the lines
  TVector3 CalculateLineMinDistanceMidpoint(
      const TVector3& pos1, const TVector3& dir1,
      const TVector3& pos2, const TVector3& dir2);

  // Helper function to calculate Pandora vertex position from two particles
  // Uses line intersection algorithm to find closest approach between two particle tracks
  // Returns: 1 if used simple average (lines parallel), 0 if used line intersection
  int CalculatePandoraVertexPosition(AnaParticlePD* particle1, AnaParticlePD* particle2, Float_t* position);

}

#endif