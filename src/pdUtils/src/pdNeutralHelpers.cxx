#include "pdNeutralHelpers.hxx"
#include <cmath>

namespace pdNeutralHelpers {

//***************************************************************
// Helper function to calculate Pandora vertex position from two particles
// Returns: 1 if used simple average (else branch), 0 if used line intersection
//***************************************************************
int CalculatePandoraVertexPosition(AnaParticlePD* particle1, AnaParticlePD* particle2, Float_t* position) {
    if (!particle1 || !particle2) {
        for (int i = 0; i < 3; i++) position[i] = -999.;
        return -999; // Default to "just average" for invalid inputs
    }

    TVector3 pos1(particle1->PositionStart[0], particle1->PositionStart[1], particle1->PositionStart[2]);
    TVector3 dir1(particle1->DirectionStart[0], particle1->DirectionStart[1], particle1->DirectionStart[2]);
    TVector3 pos2(particle2->PositionStart[0], particle2->PositionStart[1], particle2->PositionStart[2]);
    TVector3 dir2(particle2->DirectionStart[0], particle2->DirectionStart[1], particle2->DirectionStart[2]);

    // Normalize directions
    if (dir1.Mag() > 0) dir1 = dir1.Unit();
    if (dir2.Mag() > 0) dir2 = dir2.Unit();

    // Find closest point between two lines
    TVector3 w0 = pos1 - pos2;
    double a = dir1.Dot(dir1);
    double b = dir1.Dot(dir2);
    double c = dir2.Dot(dir2);
    double d = dir1.Dot(w0);
    double e = dir2.Dot(w0);
    double denom = a*c - b*b;

    TVector3 pandoraVertex;
    int isJustAverage;
    if (fabs(denom) > 1e-6) {
        double s = (b*e - c*d) / denom;
        double t = (a*e - b*d) / denom;
        TVector3 p1 = pos1 + s * dir1;
        TVector3 p2 = pos2 + t * dir2;
        pandoraVertex = 0.5 * (p1 + p2);
        isJustAverage = 0; // Used line intersection
    } else {
        pandoraVertex = 0.5 * (pos1 + pos2);
        isJustAverage = 1; // Used simple average (else branch)
    }

    position[0] = pandoraVertex.X();
    position[1] = pandoraVertex.Y();
    position[2] = pandoraVertex.Z();

    return isJustAverage;
}

//***************************************************************
TVector3 CalculateLineMinDistanceMidpoint(
    const TVector3& pos1, const TVector3& dir1,
    const TVector3& pos2, const TVector3& dir2) {
//***************************************************************

  // Calculate minimum distance point between two skew lines
  // Line 1: P1(s) = pos1 + s * dir1
  // Line 2: P2(t) = pos2 + t * dir2

  // Vector connecting the two line origins
  TVector3 w0 = pos1 - pos2;

  // Dot products needed for the formula
  double a = dir1.Dot(dir1);  // |dir1|^2
  double b = dir1.Dot(dir2);  // dir1 · dir2
  double c = dir2.Dot(dir2);  // |dir2|^2
  double d = dir1.Dot(w0);    // dir1 · w0
  double e = dir2.Dot(w0);    // dir2 · w0

  double denom = a * c - b * b;

  // Default to midpoint of line origins if lines are parallel
  if (fabs(denom) < 1e-10) {
    return 0.5 * (pos1 + pos2);
  }

  // Calculate parameters for closest points
  double sc = (b * e - c * d) / denom;
  double tc = (a * e - b * d) / denom;

  // Get the two closest points
  TVector3 P1 = pos1 + sc * dir1;
  TVector3 P2 = pos2 + tc * dir2;

  // Return midpoint
  return 0.5 * (P1 + P2);
}

} // namespace pdNeutralHelpers

