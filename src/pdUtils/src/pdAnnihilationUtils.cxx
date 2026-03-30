#include "pdAnnihilationUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <algorithm>
#include <unordered_set>

namespace pdAnnihilationUtils {

//***************************************************************
void FillPositionPandora(AnaAnnihilationVertexPD* vertex) {
//***************************************************************
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* p1 = vertex->Particles[0];
  AnaParticlePD* p2 = vertex->Particles[1];
  if (!p1 || !p2) return;

  TVector3 x1(p1->PositionStart[0], p1->PositionStart[1], p1->PositionStart[2]);
  TVector3 x2(p2->PositionStart[0], p2->PositionStart[1], p2->PositionStart[2]);
  TVector3 d1(p1->DirectionStart[0], p1->DirectionStart[1], p1->DirectionStart[2]);
  TVector3 d2(p2->DirectionStart[0], p2->DirectionStart[1], p2->DirectionStart[2]);

  const double eps = 1e-10;
  if (d1.Mag2() < eps || d2.Mag2() < eps) {
    TVector3 fallback = 0.5 * (x1 + x2);
    vertex->PositionPandora[0] = fallback.X();
    vertex->PositionPandora[1] = fallback.Y();
    vertex->PositionPandora[2] = fallback.Z();
    vertex->ClosestPointPandora1[0] = x1.X();
    vertex->ClosestPointPandora1[1] = x1.Y();
    vertex->ClosestPointPandora1[2] = x1.Z();
    vertex->ClosestPointPandora2[0] = x2.X();
    vertex->ClosestPointPandora2[1] = x2.Y();
    vertex->ClosestPointPandora2[2] = x2.Z();
    vertex->MinimumDistancePandora = static_cast<float>((x1 - x2).Mag());
    return;
  }

  d1 = d1.Unit();
  d2 = d2.Unit();
  TVector3 w0 = x1 - x2;

  const double a = d1.Dot(d1);
  const double b = d1.Dot(d2);
  const double c = d2.Dot(d2);
  const double d = d1.Dot(w0);
  const double e = d2.Dot(w0);
  const double den = a * c - b * b;

  if (den < eps && den > -eps) {
    TVector3 fallback = 0.5 * (x1 + x2);
    vertex->PositionPandora[0] = fallback.X();
    vertex->PositionPandora[1] = fallback.Y();
    vertex->PositionPandora[2] = fallback.Z();
    vertex->ClosestPointPandora1[0] = x1.X();
    vertex->ClosestPointPandora1[1] = x1.Y();
    vertex->ClosestPointPandora1[2] = x1.Z();
    vertex->ClosestPointPandora2[0] = x2.X();
    vertex->ClosestPointPandora2[1] = x2.Y();
    vertex->ClosestPointPandora2[2] = x2.Z();
    vertex->MinimumDistancePandora = static_cast<float>((x1 - x2).Mag());
    return;
  }

  const double t1 = (b * e - c * d) / den;
  const double t2 = (a * e - b * d) / den;
  TVector3 c1 = x1 + t1 * d1;
  TVector3 c2 = x2 + t2 * d2;
  TVector3 intersection = 0.5 * (c1 + c2);
  const float minDistance = static_cast<float>((c1 - c2).Mag());

  vertex->PositionPandora[0] = intersection.X();
  vertex->PositionPandora[1] = intersection.Y();
  vertex->PositionPandora[2] = intersection.Z();
  vertex->ClosestPointPandora1[0] = c1.X();
  vertex->ClosestPointPandora1[1] = c1.Y();
  vertex->ClosestPointPandora1[2] = c1.Z();
  vertex->ClosestPointPandora2[0] = c2.X();
  vertex->ClosestPointPandora2[1] = c2.Y();
  vertex->ClosestPointPandora2[2] = c2.Z();
  vertex->MinimumDistancePandora = minDistance;
}

//***************************************************************
void FillPositionFit(AnaAnnihilationVertexPD* vertex, double trackFitLength, double trackFitDistanceFromStart) {
//***************************************************************
  if (!vertex || vertex->Particles.size() < 2) return;

  AnaParticlePD* p1 = vertex->Particles[0];
  AnaParticlePD* p2 = vertex->Particles[1];
  if (!p1 || !p2) return;

  TVector3 s1(p1->PositionStart[0], p1->PositionStart[1], p1->PositionStart[2]);
  TVector3 s2(p2->PositionStart[0], p2->PositionStart[1], p2->PositionStart[2]);

  std::vector<double> fit1;
  std::vector<double> fit2;
  pdAnaUtils::ExtrapolateTrack(p1, fit1, trackFitLength, true, trackFitDistanceFromStart);
  pdAnaUtils::ExtrapolateTrack(p2, fit2, trackFitLength, true, trackFitDistanceFromStart);

  TVector3 d1(p1->DirectionStart[0], p1->DirectionStart[1], p1->DirectionStart[2]);
  TVector3 d2(p2->DirectionStart[0], p2->DirectionStart[1], p2->DirectionStart[2]);
  if (fit1.size() >= 6) d1.SetXYZ(fit1[3], fit1[4], fit1[5]);
  if (fit2.size() >= 6) d2.SetXYZ(fit2[3], fit2[4], fit2[5]);

  const double eps = 1e-10;
  if (d1.Mag2() < eps || d2.Mag2() < eps) {
    TVector3 fallback = 0.5 * (s1 + s2);
    vertex->PositionFit[0] = fallback.X();
    vertex->PositionFit[1] = fallback.Y();
    vertex->PositionFit[2] = fallback.Z();
    vertex->ClosestPointFit1[0] = s1.X();
    vertex->ClosestPointFit1[1] = s1.Y();
    vertex->ClosestPointFit1[2] = s1.Z();
    vertex->ClosestPointFit2[0] = s2.X();
    vertex->ClosestPointFit2[1] = s2.Y();
    vertex->ClosestPointFit2[2] = s2.Z();
    vertex->MinimumDistanceFit = static_cast<float>((s1 - s2).Mag());
    return;
  }
  d1 = d1.Unit();
  d2 = d2.Unit();

  std::vector<double> line1(6), line2(6);
  line1[0] = s1.X();     line1[1] = s1.Y();     line1[2] = s1.Z();
  line1[3] = d1.X();     line1[4] = d1.Y();     line1[5] = d1.Z();
  line2[0] = s2.X();     line2[1] = s2.Y();     line2[2] = s2.Z();
  line2[3] = d2.X();     line2[4] = d2.Y();     line2[5] = d2.Z();

  TVector3 closest1, closest2;
  const float minDistance = static_cast<float>(pdAnaUtils::FindClosestPointsBetweenLines(line1, line2, closest1, closest2));
  TVector3 fitPos = 0.5 * (closest1 + closest2);

  vertex->PositionFit[0] = fitPos.X();
  vertex->PositionFit[1] = fitPos.Y();
  vertex->PositionFit[2] = fitPos.Z();
  vertex->ClosestPointFit1[0] = closest1.X();
  vertex->ClosestPointFit1[1] = closest1.Y();
  vertex->ClosestPointFit1[2] = closest1.Z();
  vertex->ClosestPointFit2[0] = closest2.X();
  vertex->ClosestPointFit2[1] = closest2.Y();
  vertex->ClosestPointFit2[2] = closest2.Z();
  vertex->MinimumDistanceFit = minDistance;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> FilterVerticesByMinimumDistanceFit(std::vector<AnaAnnihilationVertexPD*>& vertices) {
//***************************************************************
  std::sort(vertices.begin(), vertices.end(),
            [](const AnaAnnihilationVertexPD* a, const AnaAnnihilationVertexPD* b) {
              return a->MinimumDistanceFit < b->MinimumDistanceFit;
            });

  std::unordered_set<AnaParticlePD*> usedParticles;
  std::vector<AnaAnnihilationVertexPD*> filtered;
  filtered.reserve(vertices.size());

  for (AnaAnnihilationVertexPD* vtx : vertices) {
    bool overlaps = false;
    for (AnaParticlePD* p : vtx->Particles) {
      if (usedParticles.find(p) != usedParticles.end()) {
        overlaps = true;
        break;
      }
    }

    if (overlaps) {
      delete vtx;
      continue;
    }

    filtered.push_back(vtx);
    for (AnaParticlePD* p : vtx->Particles) {
      usedParticles.insert(p);
    }
  }

  return filtered;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxDaughterDistance) {
//***************************************************************
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart =
      ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
  const int minCollectionHitsPerDaughter =
      ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexMinCollectionHits");

  std::vector<AnaAnnihilationVertexPD*> reconstructedVertices;
  int vertexID = 0;

  for (int i = 0; i < nParts; ++i) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if (!daughter1) continue;
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) continue;
    if (daughter1->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;

    for (int j = i + 1; j < nParts; ++j) {
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if (!daughter2) continue;
      if (daughter1->ParentID != daughter2->ParentID) continue;
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) continue;
      if (daughter2->NHitsPerPlane[2] <= minCollectionHitsPerDaughter) continue;

      TVector3 s1(daughter1->PositionStart[0], daughter1->PositionStart[1], daughter1->PositionStart[2]);
      TVector3 s2(daughter2->PositionStart[0], daughter2->PositionStart[1], daughter2->PositionStart[2]);
      float distance = (s1 - s2).Mag();
      if (distance > maxDaughterDistance) continue;

      AnaAnnihilationVertexPD* reconstructedVertex = new AnaAnnihilationVertexPD();
      reconstructedVertex->OriginalDistance = distance;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);
      reconstructedVertex->NParticles = reconstructedVertex->Particles.size();
      FillPositionPandora(reconstructedVertex);
      FillPositionFit(reconstructedVertex, trackFitLength, trackFitDistanceFromStart);
      reconstructedVertices.push_back(reconstructedVertex);
    }
  }

  return FilterVerticesByMinimumDistanceFit(reconstructedVertices);
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance) {
//***************************************************************
  return CreateVerticesCommon(event, maxDaughterDistance);
}

} // namespace pdAnnihilationUtils
