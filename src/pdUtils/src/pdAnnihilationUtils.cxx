#include "pdAnnihilationUtils.hxx"
#include "pdNeutralHelpers.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdKalman.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"
#include "TMinuit.h"
#include "TVector3.h"
#include "TMath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>

// Namespace for internal structures and functions
namespace {
  // Structure to hold hit data for TMinuit FCN
  struct VertexFitData {
    std::vector<TVector3> hits1;  // Hits from first track
    std::vector<TVector3> hits2;  // Hits from second track
  };

  // Global pointer to fit data (needed for TMinuit callback)
  VertexFitData* gFitData = nullptr;

  constexpr Int_t kTrajectoryDirHistBins = 60;
  constexpr Float_t kTrajectoryDirHistMin = -1.f;
  constexpr Float_t kTrajectoryDirHistMax = 1.f;

  struct TrajectoryDirectionResult {
    TVector3 mpv;
    Int_t nPoints = 0;

    TrajectoryDirectionResult() :
      mpv(-999.f, -999.f, -999.f) {}
  };

  TrajectoryDirectionResult ComputeTrajectoryDirectionMPV(const std::vector<AnaTrajectoryPointPD>& points) {
    TrajectoryDirectionResult result;
    if (points.empty()) {
      return result;
    }

    const Float_t binWidth = (kTrajectoryDirHistMax - kTrajectoryDirHistMin) / kTrajectoryDirHistBins;
    std::array<Int_t, kTrajectoryDirHistBins> countsX{};
    std::array<Int_t, kTrajectoryDirHistBins> countsY{};
    std::array<Int_t, kTrajectoryDirHistBins> countsZ{};

    auto fillBin = [&](Float_t value, std::array<Int_t, kTrajectoryDirHistBins>& counts) {
      if (!std::isfinite(value)) return;
      Float_t clamped = std::max(kTrajectoryDirHistMin, std::min(kTrajectoryDirHistMax - 1e-6f, value));
      Int_t bin = static_cast<Int_t>((clamped - kTrajectoryDirHistMin) / binWidth);
      bin = std::max(0, std::min(bin, kTrajectoryDirHistBins - 1));
      counts[bin]++;
    };

    for (const auto& point : points) {
      TVector3 dir = point.Direction;
      if (!std::isfinite(dir.X()) || !std::isfinite(dir.Y()) || !std::isfinite(dir.Z())) {
        continue;
      }
      if (dir.Mag() <= 0) {
        continue;
      }
      dir = dir.Unit();
      ++result.nPoints;

      fillBin(dir.X(), countsX);
      fillBin(dir.Y(), countsY);
      fillBin(dir.Z(), countsZ);
    }

    if (result.nPoints <= 0) {
    result.mpv.SetXYZ(-999.f, -999.f, -999.f);
      return result;
    }

    auto binToValue = [&](const std::array<Int_t, kTrajectoryDirHistBins>& counts) -> Float_t {
      auto it = std::max_element(counts.begin(), counts.end());
      if (it == counts.end() || *it == 0) {
        return -999.f;
      }
      Int_t bin = static_cast<Int_t>(std::distance(counts.begin(), it));
      return kTrajectoryDirHistMin + (bin + 0.5f) * binWidth;
    };

    result.mpv.SetXYZ(binToValue(countsX),
                      binToValue(countsY),
                      binToValue(countsZ));
    return result;
  }

  void EnsureTrajectoryDirectionComputed(AnaParticlePD* particle) {
    if (!particle) {
      return;
    }
    if (particle->TrajectoryDirectionNPoints > 0) {
      return;
    }

    TrajectoryDirectionResult trajResult = ComputeTrajectoryDirectionMPV(particle->TrjPoints);
    particle->TrajectoryDirection = trajResult.mpv;
    particle->TrajectoryDirectionNPoints = trajResult.nPoints;
  }

  //***************************************************************
  // TMinuit FCN function for vertex fitting
  //***************************************************************
  void VertexFitFCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // Suppress compiler warnings
    (void)npar;
    (void)gin;
    (void)iflag;

    if (!gFitData) {
      f = 1e10;
      return;
    }

    // Extract parameters
    // par[0-2]: vertex position (vx, vy, vz)
    // par[3-5]: direction 1 (dx1, dy1, dz1)
    // par[6-8]: direction 2 (dx2, dy2, dz2)
    TVector3 vertex(par[0], par[1], par[2]);
    TVector3 dir1(par[3], par[4], par[5]);
    TVector3 dir2(par[6], par[7], par[8]);

    // Normalize direction vectors
    if (dir1.Mag() > 0) dir1 = dir1.Unit();
    else {
      f = 1e10;
      return;
    }

    if (dir2.Mag() > 0) dir2 = dir2.Unit();
    else {
      f = 1e10;
      return;
    }

    // Calculate sum of squared distances along track direction
    double chi2 = 0.0;

    // Process hits from track 1
    for (const auto& hit : gFitData->hits1) {
      // Vector from vertex to hit
      TVector3 toHit = hit - vertex;

      // Project onto direction to find parameter t
      double t = toHit.Dot(dir1);

      // Calculate projection point on line
      TVector3 projection = vertex + t * dir1;

      // Calculate squared distance
      double dist2 = (hit - projection).Mag2();
      chi2 += dist2;
    }

    // Process hits from track 2
    for (const auto& hit : gFitData->hits2) {
      // Vector from vertex to hit
      TVector3 toHit = hit - vertex;

      // Project onto direction to find parameter t
      double t = toHit.Dot(dir2);

      // Calculate projection point on line
      TVector3 projection = vertex + t * dir2;

      // Calculate squared distance
      double dist2 = (hit - projection).Mag2();
      chi2 += dist2;
    }

    f = chi2;
  }
} // anonymous namespace

namespace pdAnnihilationUtils {

namespace {
  // Effective score for vertex ordering: lower is better.
  // Optional terms (all add to score, so higher penalty = worse):
  // - FilterVerticesFavourTwoPions: mass near K0 + pion PID (weak discriminator; use with care).
  // - FilterVerticesFavourGeometry: opening angle (prefer small), min daughter length (prefer long),
  //   min daughter nhits (prefer many). Based on signal vs background diagnostics.
  double GetEffectiveVertexScore(const AnaAnnihilationVertexPD* v) {
    double base = (v->Score > -900 && v->Score < 1e6) ? v->Score : 1e6;
    double total = base;

    // --- Geometry-based terms (diagnostics: signal has smaller opening, longer daughters, more hits) ---
    if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesFavourGeometry") &&
        ND::params().GetParameterI("neutralKaonAnalysis.FilterVerticesFavourGeometry") == 1) {
      double openWeight = 0.01;  // per degree
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesOpeningAngleWeight")) {
        openWeight = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesOpeningAngleWeight");
      }
      double lengthWeight = 0.02;  // per cm below threshold
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesMinDaughterLengthWeight")) {
        lengthWeight = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesMinDaughterLengthWeight");
      }
      double lengthThreshold = 40.0;  // cm
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesMinDaughterLengthThreshold")) {
        lengthThreshold = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesMinDaughterLengthThreshold");
      }
      double nhitsWeight = 0.005;  // per hit below threshold
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesMinDaughterNhitsWeight")) {
        nhitsWeight = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesMinDaughterNhitsWeight");
      }
      double nhitsThreshold = 50.0;
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesMinDaughterNhitsThreshold")) {
        nhitsThreshold = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesMinDaughterNhitsThreshold");
      }

      if (v->OpeningAngle > -900 && v->OpeningAngle < 1e3) {
        total += openWeight * v->OpeningAngle;  // degrees; signal median ~21, bkg ~56
      }
      if (v->NParticles >= 2) {
        AnaParticlePD* p1 = static_cast<AnaParticlePD*>(v->Particles[0]);
        AnaParticlePD* p2 = static_cast<AnaParticlePD*>(v->Particles[1]);
        if (p1 && p2) {
          double l1 = (p1->Length > -900 && p1->Length < 1e4) ? p1->Length : 0;
          double l2 = (p2->Length > -900 && p2->Length < 1e4) ? p2->Length : 0;
          double minLen = std::min(l1, l2);
          if (minLen >= 0) {
            total += lengthWeight * std::max(0.0, lengthThreshold - minLen);
          }
          int n1 = static_cast<int>(p1->Hits[2].size());
          int n2 = static_cast<int>(p2->Hits[2].size());
          int minNhits = std::min(n1, n2);
          if (minNhits >= 0) {
            total += nhitsWeight * std::max(0.0, nhitsThreshold - static_cast<double>(minNhits));
          }
        }
      }
    }

    // --- Mass + PID terms (optional; weak discriminator when background is also pi+ pi-) ---
    if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesFavourTwoPions") &&
        ND::params().GetParameterI("neutralKaonAnalysis.FilterVerticesFavourTwoPions") == 1) {
      double massWeight = 1e-5;
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesK0MassWeight")) {
        massWeight = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesK0MassWeight");
      }
      double pidWeight = 0.01;
      if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesPionChi2Weight")) {
        pidWeight = ND::params().GetParameterD("neutralKaonAnalysis.FilterVerticesPionChi2Weight");
      }
      const double M_K0_MeV = 497.611;
      double E_MeV = v->Energy;
      double px = v->Momentum[0], py = v->Momentum[1], pz = v->Momentum[2];
      double p2_MeV2 = 1e6 * (px*px + py*py + pz*pz);
      double M2 = E_MeV * E_MeV - p2_MeV2;
      double M_MeV = (M2 > 0) ? std::sqrt(M2) : 0;
      total += massWeight * (M_MeV - M_K0_MeV) * (M_MeV - M_K0_MeV);
      if (v->NParticles >= 2 && pidWeight > 0) {
        AnaParticlePD* p1 = static_cast<AnaParticlePD*>(v->Particles[0]);
        AnaParticlePD* p2 = static_cast<AnaParticlePD*>(v->Particles[1]);
        if (p1 && p2) {
          std::pair<double, int> c1 = pdAnaUtils::Chi2PID(*p1, 211);
          std::pair<double, int> c2 = pdAnaUtils::Chi2PID(*p2, 211);
          double ndf1 = (c1.second > 0) ? c1.first / c1.second : 9999.0;
          double ndf2 = (c2.second > 0) ? c2.first / c2.second : 9999.0;
          total += pidWeight * std::max(ndf1, ndf2);
        }
      }
    }

    return total;
  }
} // anonymous namespace

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> FilterVerticesByScore(std::vector<AnaAnnihilationVertexPD*>& vertices) {
    //***************************************************************

      static bool once = true;
      if (once && !vertices.empty()) {
        once = false;
        if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesFavourGeometry") &&
            ND::params().GetParameterI("neutralKaonAnalysis.FilterVerticesFavourGeometry") == 1) {
          std::cout << "[FilterVerticesByScore] FilterVerticesFavourGeometry=1: ordering by score + opening angle + min daughter length/nhits" << std::endl;
        }
        if (ND::params().HasParameter("neutralKaonAnalysis.FilterVerticesFavourTwoPions") &&
            ND::params().GetParameterI("neutralKaonAnalysis.FilterVerticesFavourTwoPions") == 1) {
          std::cout << "[FilterVerticesByScore] FilterVerticesFavourTwoPions=1: ordering also by mass near K0 + pion PID" << std::endl;
        }
      }
      // Sort vertices by effective score (ascending - lower is better).
      // Optionally favour two-pion vertices (mass near K0, pion-like PID) via neutralKaonAnalysis.FilterVerticesFavourTwoPions.
      std::sort(vertices.begin(), vertices.end(),
        [](const AnaAnnihilationVertexPD* a, const AnaAnnihilationVertexPD* b) {
          double sa = GetEffectiveVertexScore(a);
          double sb = GetEffectiveVertexScore(b);
          return sa < sb;
        });

      // Track which particles have been used
      std::set<AnaParticlePD*> usedParticles;

      // Output vector for selected vertices
      std::vector<AnaAnnihilationVertexPD*> selectedVertices;
      std::set<AnaAnnihilationVertexPD*> selectedSet;  // For quick lookup

      // Iterate through sorted vertices
      for (AnaAnnihilationVertexPD* vertex : vertices) {
        if (!vertex || vertex->NParticles < 2) continue;

        AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
        AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

        if (!part1 || !part2) continue;

        // Check if either particle is already used
        if (usedParticles.find(part1) == usedParticles.end() &&
            usedParticles.find(part2) == usedParticles.end()) {
          // Neither particle is used, so we can keep this vertex
          selectedVertices.push_back(vertex);
          selectedSet.insert(vertex);

          // Mark particles as used
          usedParticles.insert(part1);
          usedParticles.insert(part2);
        }
      }

      // Delete vertices that were not selected
      for (AnaAnnihilationVertexPD* vertex : vertices) {
        if (selectedSet.find(vertex) == selectedSet.end()) {
          delete vertex;
        }
      }

      return selectedVertices;
    }

//***************************************************************
double FindVertexPositionWithFit(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!part1 || !part2) {
    return -999.0;
  }

  // Get reference positions for both particles
  TVector3 pos1 = pdAnaUtils::DefinePosition(part1, true);  // Use start position
  TVector3 pos2 = pdAnaUtils::DefinePosition(part2, true);  // Use start position

  // Check if positions are valid
  if (pos1.X() < -900 || pos2.X() < -900) {
    return -999.0;
  }

  // Collect hits within trackFitLength from each particle
  VertexFitData fitData;

  // Collect hits for particle 1
  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part1->Hits[plane].size(); i++) {
      AnaHitPD& hit = part1->Hits[plane][i];
      if (hit.Position.Z() != -999) {
        TVector3 hitPos = hit.Position;
        double distance = (hitPos - pos1).Mag();
        if (distance <= trackFitLength) {
          fitData.hits1.push_back(hitPos);
        }
      }
    }
  }

  // Collect hits for particle 2
  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part2->Hits[plane].size(); i++) {
      AnaHitPD& hit = part2->Hits[plane][i];
      if (hit.Position.Z() != -999) {
        TVector3 hitPos = hit.Position;
        double distance = (hitPos - pos2).Mag();
        if (distance <= trackFitLength) {
          fitData.hits2.push_back(hitPos);
        }
      }
    }
  }

  // Need sufficient hits for fitting
  if (fitData.hits1.size() < 2 || fitData.hits2.size() < 2) {
    return -999.0;
  }

  // Set global fit data pointer
  gFitData = &fitData;

  // Initialize TMinuit with 9 parameters
  TMinuit minuit(9);
  minuit.SetPrintLevel(-1);  // Suppress output

  // Set FCN function
  minuit.SetFCN(VertexFitFCN);

  // Calculate initial vertex position as midpoint
  double vx_init = 0.5 * (pos1.X() + pos2.X());
  double vy_init = 0.5 * (pos1.Y() + pos2.Y());
  double vz_init = 0.5 * (pos1.Z() + pos2.Z());

  // Get initial directions from particles
  TVector3 dir1_init(part1->DirectionStart[0], part1->DirectionStart[1], part1->DirectionStart[2]);
  TVector3 dir2_init(part2->DirectionStart[0], part2->DirectionStart[1], part2->DirectionStart[2]);

  // Normalize initial directions
  if (dir1_init.Mag() > 0) dir1_init = dir1_init.Unit();
  if (dir2_init.Mag() > 0) dir2_init = dir2_init.Unit();

  // Set parameters
  Double_t step = 0.1;
  minuit.DefineParameter(0, "vx", vx_init, step, 0, 0);
  minuit.DefineParameter(1, "vy", vy_init, step, 0, 0);
  minuit.DefineParameter(2, "vz", vz_init, step, 0, 0);
  minuit.DefineParameter(3, "dx1", dir1_init.X(), step, 0, 0);
  minuit.DefineParameter(4, "dy1", dir1_init.Y(), step, 0, 0);
  minuit.DefineParameter(5, "dz1", dir1_init.Z(), step, 0, 0);
  minuit.DefineParameter(6, "dx2", dir2_init.X(), step, 0, 0);
  minuit.DefineParameter(7, "dy2", dir2_init.Y(), step, 0, 0);
  minuit.DefineParameter(8, "dz2", dir2_init.Z(), step, 0, 0);

  // Run MIGRAD minimization
  minuit.Migrad();

  // Get fitted parameters
  Double_t vx, vy, vz, dx1, dy1, dz1, dx2, dy2, dz2;
  Double_t err;
  minuit.GetParameter(0, vx, err);
  minuit.GetParameter(1, vy, err);
  minuit.GetParameter(2, vz, err);
  minuit.GetParameter(3, dx1, err);
  minuit.GetParameter(4, dy1, err);
  minuit.GetParameter(5, dz1, err);
  minuit.GetParameter(6, dx2, err);
  minuit.GetParameter(7, dy2, err);
  minuit.GetParameter(8, dz2, err);

  // Get chi2 value
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  // Store fitted vertex position in PositionFit
  vertex->PositionFit[0] = vx;
  vertex->PositionFit[1] = vy;
  vertex->PositionFit[2] = vz;

  // Calculate Pandora-based vertex position and set flag
  vertex->IsJustAverage = pdNeutralHelpers::CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

  // Copy Pandora position to main Position
  vertex->Position[0] = vertex->PositionPandora[0];
  vertex->Position[1] = vertex->PositionPandora[1];
  vertex->Position[2] = vertex->PositionPandora[2];

  // Store Pandora-based closest points for downstream visualization
  if (part1 && part2) {
    TVector3 pandoraPos1(part1->PositionStart[0], part1->PositionStart[1], part1->PositionStart[2]);
    TVector3 pandoraDir1(part1->DirectionStart[0], part1->DirectionStart[1], part1->DirectionStart[2]);
    TVector3 pandoraPos2(part2->PositionStart[0], part2->PositionStart[1], part2->PositionStart[2]);
    TVector3 pandoraDir2(part2->DirectionStart[0], part2->DirectionStart[1], part2->DirectionStart[2]);

    if (pandoraDir1.Mag() > 0) pandoraDir1 = pandoraDir1.Unit();
    if (pandoraDir2.Mag() > 0) pandoraDir2 = pandoraDir2.Unit();

    TVector3 w0 = pandoraPos1 - pandoraPos2;
    double a = pandoraDir1.Dot(pandoraDir1);
    double b = pandoraDir1.Dot(pandoraDir2);
    double c = pandoraDir2.Dot(pandoraDir2);
    double d = pandoraDir1.Dot(w0);
    double e = pandoraDir2.Dot(w0);
    double denom = a * c - b * b;

    if (fabs(denom) > 1e-6) {
      double s = (b * e - c * d) / denom;
      double t = (a * e - b * d) / denom;

      TVector3 closest1 = pandoraPos1 + s * pandoraDir1;
      TVector3 closest2 = pandoraPos2 + t * pandoraDir2;

      vertex->ClosestPoint1[0] = closest1.X();
      vertex->ClosestPoint1[1] = closest1.Y();
      vertex->ClosestPoint1[2] = closest1.Z();
      vertex->ClosestPoint2[0] = closest2.X();
      vertex->ClosestPoint2[1] = closest2.Y();
      vertex->ClosestPoint2[2] = closest2.Z();
    } else {
      vertex->ClosestPoint1[0] = vertex->Position[0];
      vertex->ClosestPoint1[1] = vertex->Position[1];
      vertex->ClosestPoint1[2] = vertex->Position[2];
      vertex->ClosestPoint2[0] = vertex->Position[0];
      vertex->ClosestPoint2[1] = vertex->Position[1];
      vertex->ClosestPoint2[2] = vertex->Position[2];
    }
  } else {
    vertex->ClosestPoint1[0] = vertex->Position[0];
    vertex->ClosestPoint1[1] = vertex->Position[1];
    vertex->ClosestPoint1[2] = vertex->Position[2];
    vertex->ClosestPoint2[0] = vertex->Position[0];
    vertex->ClosestPoint2[1] = vertex->Position[1];
    vertex->ClosestPoint2[2] = vertex->Position[2];
  }

  // Calculate chi2/ndf
  int ndf = (fitData.hits1.size() + fitData.hits2.size()) * 3 - 9;  // 3 coords per hit, 9 parameters
  double chi2ndf = (ndf > 0) ? amin / ndf : amin;

  // Store score (chi2/ndf from minimization)
  vertex->Score = chi2ndf;

  // Calculate minimum distance between fitted lines
  TVector3 vertexPos(vx, vy, vz);
  TVector3 fittedDir1(dx1, dy1, dz1);
  TVector3 fittedDir2(dx2, dy2, dz2);

  // Normalize directions
  if (fittedDir1.Mag() > 0) fittedDir1 = fittedDir1.Unit();
  if (fittedDir2.Mag() > 0) fittedDir2 = fittedDir2.Unit();

  // Both lines emanate from same vertex point, so minimum distance is 0.0
  vertex->MinimumDistance = 0.0;

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = fittedDir1 + fittedDir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->DirectionFit[0] = fitDirSum.X();
  vertex->DirectionFit[1] = fitDirSum.Y();
  vertex->DirectionFit[2] = fitDirSum.Z();

  // Clean up global pointer
  gFitData = nullptr;

  return chi2ndf;
}

//***************************************************************
double FindVertexPositionGeometric(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  (void)trackFitLength; // Not used in geometric method

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // ========== SECTION 1: Fit calculations ==========
  // Call FindVertexPositionFit to get Fit results
  // This automatically fills: PositionFit, MinimumDistanceFit, ClosestPoint1Fit, ClosestPoint2Fit, FittedLineParams
  FindVertexPositionFit(vertex);

  // Calculate fitted direction from the stored line parameters
  // DirectionFit = sum of directions from the two fitted lines
  if(vertex->FittedLineParams.size() >= 2){
    // Extract line parameters for daughter 1
    std::vector<double> line1 = vertex->FittedLineParams[0];
    TVector3 dir1(line1[3], line1[4], line1[5]);
    if(dir1.Mag() > 0) dir1 = dir1.Unit();

    // Extract line parameters for daughter 2
    std::vector<double> line2 = vertex->FittedLineParams[1];
    TVector3 dir2(line2[3], line2[4], line2[5]);
    if(dir2.Mag() > 0) dir2 = dir2.Unit();

    // Calculate vertex fit direction (sum of daughter fit directions, normalized)
    TVector3 fitDirSum = dir1 + dir2;
    if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
    vertex->DirectionFit[0] = fitDirSum.X();
    vertex->DirectionFit[1] = fitDirSum.Y();
    vertex->DirectionFit[2] = fitDirSum.Z();
  } else {
    // Fitted line parameters not available, set to invalid
    vertex->DirectionFit[0] = -999.0;
    vertex->DirectionFit[1] = -999.0;
    vertex->DirectionFit[2] = -999.0;
  }

  // ========== SECTION 2: Pandora calculations ==========
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  // Calculate Pandora-based vertex position and set flag
  vertex->IsJustAverage = pdNeutralHelpers::CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

  // Calculate DirectionPandora as normalized sum of Pandora directions
  // DirectionPandora = sum of directions from the two particles' DirectionStart (Pandora directions)
  if (part1 && part2) {
    TVector3 pandoraDir1(part1->DirectionStart[0], part1->DirectionStart[1], part1->DirectionStart[2]);
    TVector3 pandoraDir2(part2->DirectionStart[0], part2->DirectionStart[1], part2->DirectionStart[2]);

    if (pandoraDir1.Mag() > 0) pandoraDir1 = pandoraDir1.Unit();
    if (pandoraDir2.Mag() > 0) pandoraDir2 = pandoraDir2.Unit();

    TVector3 pandoraDirSum = pandoraDir1 + pandoraDir2;
    if (pandoraDirSum.Mag() > 0) pandoraDirSum = pandoraDirSum.Unit();
    vertex->DirectionPandora[0] = pandoraDirSum.X();
    vertex->DirectionPandora[1] = pandoraDirSum.Y();
    vertex->DirectionPandora[2] = pandoraDirSum.Z();

    // Calculate closest points between Pandora lines
    TVector3 pandoraPos1(part1->PositionStart[0], part1->PositionStart[1], part1->PositionStart[2]);
    TVector3 pandoraPos2(part2->PositionStart[0], part2->PositionStart[1], part2->PositionStart[2]);

    TVector3 w0 = pandoraPos1 - pandoraPos2;
    double a = pandoraDir1.Dot(pandoraDir1);
    double b = pandoraDir1.Dot(pandoraDir2);
    double c = pandoraDir2.Dot(pandoraDir2);
    double d = pandoraDir1.Dot(w0);
    double e = pandoraDir2.Dot(w0);
    double denom = a * c - b * b;

    if (fabs(denom) > 1e-6) {
      double s = (b * e - c * d) / denom;
      double t = (a * e - b * d) / denom;

      TVector3 closest1 = pandoraPos1 + s * pandoraDir1;
      TVector3 closest2 = pandoraPos2 + t * pandoraDir2;

      vertex->ClosestPoint1Pandora[0] = closest1.X();
      vertex->ClosestPoint1Pandora[1] = closest1.Y();
      vertex->ClosestPoint1Pandora[2] = closest1.Z();
      vertex->ClosestPoint2Pandora[0] = closest2.X();
      vertex->ClosestPoint2Pandora[1] = closest2.Y();
      vertex->ClosestPoint2Pandora[2] = closest2.Z();

      // Calculate minimum distance between Pandora lines
      vertex->MinimumDistancePandora = static_cast<Float_t>((closest1 - closest2).Mag());
    } else {
      // Lines are parallel, use Pandora position
      vertex->ClosestPoint1Pandora[0] = vertex->PositionPandora[0];
      vertex->ClosestPoint1Pandora[1] = vertex->PositionPandora[1];
      vertex->ClosestPoint1Pandora[2] = vertex->PositionPandora[2];
      vertex->ClosestPoint2Pandora[0] = vertex->PositionPandora[0];
      vertex->ClosestPoint2Pandora[1] = vertex->PositionPandora[1];
      vertex->ClosestPoint2Pandora[2] = vertex->PositionPandora[2];
      vertex->MinimumDistancePandora = 0.0;
    }
  } else {
    // Invalid particles, use Pandora position as fallback
    vertex->ClosestPoint1Pandora[0] = vertex->PositionPandora[0];
    vertex->ClosestPoint1Pandora[1] = vertex->PositionPandora[1];
    vertex->ClosestPoint1Pandora[2] = vertex->PositionPandora[2];
    vertex->ClosestPoint2Pandora[0] = vertex->PositionPandora[0];
    vertex->ClosestPoint2Pandora[1] = vertex->PositionPandora[1];
    vertex->ClosestPoint2Pandora[2] = vertex->PositionPandora[2];
    vertex->MinimumDistancePandora = -999.0;
    vertex->DirectionPandora[0] = -999.0;
    vertex->DirectionPandora[1] = -999.0;
    vertex->DirectionPandora[2] = -999.0;
  }

  // ========== SECTION 3: Parameter-based selection ==========
  // Read parameter to decide which version to use for general variables
  int usePandora = ND::params().GetParameterI("neutralKaonAnalysis.UsePandora");

  if (usePandora == 1) {
    // Use Pandora versions for general variables
    vertex->Position[0] = vertex->PositionPandora[0];
    vertex->Position[1] = vertex->PositionPandora[1];
    vertex->Position[2] = vertex->PositionPandora[2];
    vertex->Direction[0] = vertex->DirectionPandora[0];
    vertex->Direction[1] = vertex->DirectionPandora[1];
    vertex->Direction[2] = vertex->DirectionPandora[2];
    vertex->ClosestPoint1[0] = vertex->ClosestPoint1Pandora[0];
    vertex->ClosestPoint1[1] = vertex->ClosestPoint1Pandora[1];
    vertex->ClosestPoint1[2] = vertex->ClosestPoint1Pandora[2];
    vertex->ClosestPoint2[0] = vertex->ClosestPoint2Pandora[0];
    vertex->ClosestPoint2[1] = vertex->ClosestPoint2Pandora[1];
    vertex->ClosestPoint2[2] = vertex->ClosestPoint2Pandora[2];
    vertex->MinimumDistance = vertex->MinimumDistancePandora;
  } else {
    // Use Fit versions for general variables
    vertex->Position[0] = vertex->PositionFit[0];
    vertex->Position[1] = vertex->PositionFit[1];
    vertex->Position[2] = vertex->PositionFit[2];
    vertex->Direction[0] = vertex->DirectionFit[0];
    vertex->Direction[1] = vertex->DirectionFit[1];
    vertex->Direction[2] = vertex->DirectionFit[2];
    vertex->ClosestPoint1[0] = vertex->ClosestPoint1Fit[0];
    vertex->ClosestPoint1[1] = vertex->ClosestPoint1Fit[1];
    vertex->ClosestPoint1[2] = vertex->ClosestPoint1Fit[2];
    vertex->ClosestPoint2[0] = vertex->ClosestPoint2Fit[0];
    vertex->ClosestPoint2[1] = vertex->ClosestPoint2Fit[1];
    vertex->ClosestPoint2[2] = vertex->ClosestPoint2Fit[2];
    vertex->MinimumDistance = vertex->MinimumDistanceFit;
  }

  // ========== SECTION 4: Final assignment ==========
  // Always use Pandora minimum distance for Score when valid, so that FilterVerticesByScore
  // keeps the same set of vertices (and thus the same truth matching) regardless of UsePandora.
  // UsePandora only chooses which Position/Direction are written to the tree for resolution studies.
  vertex->Score = (vertex->MinimumDistancePandora > -900.f) ? vertex->MinimumDistancePandora : vertex->MinimumDistance;

  return vertex->MinimumDistance;
}

//***************************************************************
double FindVertexPositionKalman(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!daughter1 || !daughter2) {
    return -999.0;
  }

  // Check if daughter1 and daughter2 are close enough
  TVector3 pos1 = pdAnaUtils::DefinePosition(daughter1);
  TVector3 pos2 = pdAnaUtils::DefinePosition(daughter2);

  // Check if positions are valid
  if (pos1.X() < -900 || pos2.X() < -900) {
    return -999.0;
  }

  // Get seed position from geometric approach (midpoint)
  TVector3 seedPosition = 0.5 * (pos1 + pos2);

  // Convert particles to track states
  pdKalman::TrackState track1 = pdKalman::ParticleToTrackState(daughter1, trackFitLength, true);
  pdKalman::TrackState track2 = pdKalman::ParticleToTrackState(daughter2, trackFitLength, true);

  // Perform Kalman vertex fit
  pdKalman::VertexState vertexState = pdKalman::FitVertex(track1, track2, seedPosition);

  // Check if fit quality is reasonable
  double maxKalmanChi2Ndf = ND::params().GetParameterD("neutralKaonAnalysis.MaxKalmanChi2Ndf");
  double chi2ndf = (vertexState.ndf > 0) ? vertexState.chi2 / vertexState.ndf : 9999.0;
  if (chi2ndf > maxKalmanChi2Ndf) {
    // Mark vertex as invalid by setting fitted position to -999
    vertex->PositionFit[0] = -999.0;
    vertex->PositionFit[1] = -999.0;
    vertex->PositionFit[2] = -999.0;
    // Also mark main Position as invalid
    vertex->Position[0] = -999.0;
    vertex->Position[1] = -999.0;
    vertex->Position[2] = -999.0;
    vertex->Score = 9999.0;
    return 9999.0;
  }

  // Set fitted vertex position in PositionFit
  vertex->PositionFit[0] = vertexState.position.X();
  vertex->PositionFit[1] = vertexState.position.Y();
  vertex->PositionFit[2] = vertexState.position.Z();

  // Calculate Pandora-based vertex position and set flag
  vertex->IsJustAverage = pdNeutralHelpers::CalculatePandoraVertexPosition(daughter1, daughter2, vertex->PositionPandora);

  // Copy Pandora position to main Position
  vertex->Position[0] = vertex->PositionPandora[0];
  vertex->Position[1] = vertex->PositionPandora[1];
  vertex->Position[2] = vertex->PositionPandora[2];

  if (daughter1 && daughter2) {
    TVector3 pandoraPos1(daughter1->PositionStart[0], daughter1->PositionStart[1], daughter1->PositionStart[2]);
    TVector3 pandoraDir1(daughter1->DirectionStart[0], daughter1->DirectionStart[1], daughter1->DirectionStart[2]);
    TVector3 pandoraPos2(daughter2->PositionStart[0], daughter2->PositionStart[1], daughter2->PositionStart[2]);
    TVector3 pandoraDir2(daughter2->DirectionStart[0], daughter2->DirectionStart[1], daughter2->DirectionStart[2]);

    if (pandoraDir1.Mag() > 0) pandoraDir1 = pandoraDir1.Unit();
    if (pandoraDir2.Mag() > 0) pandoraDir2 = pandoraDir2.Unit();

    TVector3 w0 = pandoraPos1 - pandoraPos2;
    double a = pandoraDir1.Dot(pandoraDir1);
    double b = pandoraDir1.Dot(pandoraDir2);
    double c = pandoraDir2.Dot(pandoraDir2);
    double d = pandoraDir1.Dot(w0);
    double e = pandoraDir2.Dot(w0);
    double denom = a * c - b * b;

    if (fabs(denom) > 1e-6) {
      double s = (b * e - c * d) / denom;
      double t = (a * e - b * d) / denom;

      TVector3 closest1 = pandoraPos1 + s * pandoraDir1;
      TVector3 closest2 = pandoraPos2 + t * pandoraDir2;

      vertex->ClosestPoint1[0] = closest1.X();
      vertex->ClosestPoint1[1] = closest1.Y();
      vertex->ClosestPoint1[2] = closest1.Z();
      vertex->ClosestPoint2[0] = closest2.X();
      vertex->ClosestPoint2[1] = closest2.Y();
      vertex->ClosestPoint2[2] = closest2.Z();
    } else {
      vertex->ClosestPoint1[0] = vertex->Position[0];
      vertex->ClosestPoint1[1] = vertex->Position[1];
      vertex->ClosestPoint1[2] = vertex->Position[2];
      vertex->ClosestPoint2[0] = vertex->Position[0];
      vertex->ClosestPoint2[1] = vertex->Position[1];
      vertex->ClosestPoint2[2] = vertex->Position[2];
    }
  } else {
    vertex->ClosestPoint1[0] = vertex->Position[0];
    vertex->ClosestPoint1[1] = vertex->Position[1];
    vertex->ClosestPoint1[2] = vertex->Position[2];
    vertex->ClosestPoint2[0] = vertex->Position[0];
    vertex->ClosestPoint2[1] = vertex->Position[1];
    vertex->ClosestPoint2[2] = vertex->Position[2];
  }

  // Set score to chi2 (lower is better)
  vertex->Score = vertexState.chi2;

  // Calculate minimum distance between fitted tracks at vertex
  TVector3 vtxPos = vertexState.position;
  TVector3 track1ToVtx = vtxPos - track1.position;
  TVector3 track2ToVtx = vtxPos - track2.position;
  double dist1 = (track1ToVtx - track1.direction * track1ToVtx.Dot(track1.direction)).Mag();
  double dist2 = (track2ToVtx - track2.direction * track2ToVtx.Dot(track2.direction)).Mag();
  vertex->MinimumDistance = sqrt(dist1*dist1 + dist2*dist2);

  // Calculate vertex fit direction from Kalman track states
  TVector3 dir1 = track1.direction.Unit();
  TVector3 dir2 = track2.direction.Unit();

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = dir1 + dir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->DirectionFit[0] = fitDirSum.X();
  vertex->DirectionFit[1] = fitDirSum.Y();
  vertex->DirectionFit[2] = fitDirSum.Z();


  return vertexState.chi2;
}

//***************************************************************
bool ValidateVertex(AnaVertexPD* vertex) {
//***************************************************************

  if (!vertex) return false;

  // Check valid position (> -900 for all coordinates)
  if (vertex->Position[0] < -900 || vertex->Position[1] < -900 || vertex->Position[2] < -900) {
    return false;
  }

  // Check valid score
  if (vertex->Score < -900 || vertex->Score > 1e6) {
    return false;
  }

  // Require the two track lines to actually meet: closest-approach distance must be below threshold.
  // (PositionStart of each track can be 20+ cm apart in LArTPC because they are "first hit", not the
  // true decay point; the geometric closest approach is the physically meaningful quantity.)
  if (ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexMaxClosestApproach")) {
    double maxClosest = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexMaxClosestApproach");
    if (vertex->MinimumDistance > -900 && vertex->MinimumDistance > maxClosest) {
      return false;
    }
  }

  return true;
}

//***************************************************************
AnaTrueEquivalentVertexPD* FillTrueEquivalentVertex(AnaVertexPD* vertex) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return nullptr;
  }

  AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!daughter1 || !daughter2) {
    return nullptr;
  }

  // Create true equivalent vertex
  AnaTrueEquivalentVertexPD* trueEquivalentVertex = new AnaTrueEquivalentVertexPD();

  AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

  trueEquivalentVertex->TrueParticles.push_back(trueParticle1);
  trueEquivalentVertex->TrueParticles.push_back(trueParticle2);

  // Calculate original distance between true particles
  if (trueParticle1 && trueParticle2) {
    float trueDistance = sqrt(pow(trueParticle1->Position[0] - trueParticle2->Position[0], 2) +
                              pow(trueParticle1->Position[1] - trueParticle2->Position[1], 2) +
                              pow(trueParticle1->Position[2] - trueParticle2->Position[2], 2));
    trueEquivalentVertex->OriginalDistance = trueDistance;
  } else {
    trueEquivalentVertex->OriginalDistance = -999.0;
  }

  // Initialize position and distance
  trueEquivalentVertex->Position[0] = -999.0;
  trueEquivalentVertex->Position[1] = -999.0;
  trueEquivalentVertex->Position[2] = -999.0;
  trueEquivalentVertex->MinimumDistance = -999.0;

  // Calculate true vertex position
  FindVertexPositionFit(trueEquivalentVertex);

  // Copy PositionFit to Position (for true equivalent vertices, we use Fit results directly)
  trueEquivalentVertex->Position[0] = trueEquivalentVertex->PositionFit[0];
  trueEquivalentVertex->Position[1] = trueEquivalentVertex->PositionFit[1];
  trueEquivalentVertex->Position[2] = trueEquivalentVertex->PositionFit[2];

  // Calculate true direction (sum of daughter directions, normalized)
  if (trueParticle1 && trueParticle2) {
    TVector3 trueDir1(trueParticle1->Direction[0], trueParticle1->Direction[1], trueParticle1->Direction[2]);
    TVector3 trueDir2(trueParticle2->Direction[0], trueParticle2->Direction[1], trueParticle2->Direction[2]);
    TVector3 trueDirSum = trueDir1 + trueDir2;
    if(trueDirSum.Mag() > 0) trueDirSum = trueDirSum.Unit();
    trueEquivalentVertex->Direction[0] = trueDirSum.X();
    trueEquivalentVertex->Direction[1] = trueDirSum.Y();
    trueEquivalentVertex->Direction[2] = trueDirSum.Z();
  } else {
    trueEquivalentVertex->Direction[0] = -999.0;
    trueEquivalentVertex->Direction[1] = -999.0;
    trueEquivalentVertex->Direction[2] = -999.0;
  }

  // Copy fitted direction from reconstructed vertex
  trueEquivalentVertex->DirectionFit[0] = vertex->DirectionFit[0];
  trueEquivalentVertex->DirectionFit[1] = vertex->DirectionFit[1];
  trueEquivalentVertex->DirectionFit[2] = vertex->DirectionFit[2];

  // Initialize Pandora position
  for (int i = 0; i < 3; i++) {
    trueEquivalentVertex->PositionPandora[i] = -999;
  }

  return trueEquivalentVertex;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(
    AnaEventB& event,
    double maxDaughterDistance,
    double (*positionFinder)(AnaVertexPD*, double)) {
//***************************************************************

  // Get the array of particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // Build hash map for O(1) particle lookups by UniqueID
  std::unordered_map<Int_t, AnaParticlePD*> particleByUniqueID;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part){
      particleByUniqueID[part->UniqueID] = part;
    }
  }

  // Create reconstructed vertices
  std::vector<AnaAnnihilationVertexPD*> reconstructedVertices;
  int vertexID = 0; // Counter for unique vertex IDs

  // Set to track which particle pairs have already been used to create vertices
  std::set<std::pair<AnaParticlePD*, AnaParticlePD*>> usedPairs;

  // Get track fit length from parameters
  double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

  // Get beam particle direction for angle calculations
  TVector3 beamDir(0, 0, 1); // Default to z-axis if beam not available
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(event.Beam);
  if (beam && beam->BeamParticle) {
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    if (beamPart) {
      // Use beam particle direction (preferably end direction at TPC entrance)
      if (beamPart->DirectionEnd[0] != -999 && beamPart->DirectionEnd[1] != -999 && beamPart->DirectionEnd[2] != -999) {
        beamDir.SetXYZ(beamPart->DirectionEnd[0], beamPart->DirectionEnd[1], beamPart->DirectionEnd[2]);
      } else if (beamPart->DirectionStart[0] != -999 && beamPart->DirectionStart[1] != -999 && beamPart->DirectionStart[2] != -999) {
        beamDir.SetXYZ(beamPart->DirectionStart[0], beamPart->DirectionStart[1], beamPart->DirectionStart[2]);
      }
      beamDir = beamDir.Unit(); // Ensure normalized
    }
  }

  for(int i = 0; i < nParts; i++){
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if(!daughter1) continue;

    // Check if daughter1 has valid start position
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) {
      continue; // Skip daughter1 with invalid start positions
    }
    // Check if daughter1 has valid end position
    if (daughter1->PositionEnd[0] < -900 || daughter1->PositionEnd[1] < -900 || daughter1->PositionEnd[2] < -900) {
      continue; // Skip daughter1 with invalid end positions
    }

    // Check PID compatibility for daughter1 - skip if compatible with proton or kaon
    std::pair<double, int> chi2Proton1 = pdAnaUtils::Chi2PID(*daughter1, 2212);  // Proton PDG
    std::pair<double, int> chi2Kaon1 = pdAnaUtils::Chi2PID(*daughter1, 321);    // Kaon PDG
    double chi2ndfProton1 = (chi2Proton1.second > 0) ? chi2Proton1.first / chi2Proton1.second : 9999.0;
    double chi2ndfKaon1 = (chi2Kaon1.second > 0) ? chi2Kaon1.first / chi2Kaon1.second : 9999.0;

    if (chi2ndfProton1 < 100.0 || chi2ndfKaon1 < 20.0) {
      continue; // Skip particles compatible with proton (chi2/ndf<100) or kaon (chi2/ndf<20)
    }

    for(int j = i + 1; j < nParts; j++){
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if(!daughter2) continue;

      // Skip if daughter1 and daughter2 are the same particle
      if (daughter1 == daughter2) {
        continue;
      }

      // // Check if both particles have the same ParentID (EARLY CHECK - avoid useless computation)
      if (daughter1->ParentID != daughter2->ParentID) {
        continue; // Skip if parents are different
      }

      // Check if daughter2 has valid start position
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) {
        continue; // Skip daughter2 with invalid start positions
      }
      // Check if daughter2 has valid end position
      if (daughter2->PositionEnd[0] < -900 || daughter2->PositionEnd[1] < -900 || daughter2->PositionEnd[2] < -900) {
        continue; // Skip daughter2 with invalid end positions
      }

      // Check PID compatibility for daughter2 - skip if compatible with proton or kaon
      std::pair<double, int> chi2Proton2 = pdAnaUtils::Chi2PID(*daughter2, 2212);  // Proton PDG
      std::pair<double, int> chi2Kaon2 = pdAnaUtils::Chi2PID(*daughter2, 321);    // Kaon PDG
      double chi2ndfProton2 = (chi2Proton2.second > 0) ? chi2Proton2.first / chi2Proton2.second : 9999.0;
      double chi2ndfKaon2 = (chi2Kaon2.second > 0) ? chi2Kaon2.first / chi2Kaon2.second : 9999.0;

      if (chi2ndfProton2 < 100.0 || chi2ndfKaon2 < 20.0) {
        continue; // Skip particles compatible with proton (chi2/ndf<100) or kaon (chi2/ndf<20)
      }

      // Create a pair to check for duplicates (order-independent)
      std::pair<AnaParticlePD*, AnaParticlePD*> currentPair;
      if (daughter1 < daughter2) {
        currentPair = std::make_pair(daughter1, daughter2);
      } else {
        currentPair = std::make_pair(daughter2, daughter1);
      }

      // Check if this pair has already been used to create a vertex
      if (usedPairs.find(currentPair) != usedPairs.end()) {
        continue; // Skip if this pair has already been used
      }

      // Check if daughter1 and daughter2 are close enough using closest approach of the two
      // track lines (not just endpoints). This recovers pairs where the vertex is not at any
      // of the four endpoints but the lines still pass within maxDaughterDistance.
      TVector3 s1(daughter1->PositionStart[0], daughter1->PositionStart[1], daughter1->PositionStart[2]);
      TVector3 e1(daughter1->PositionEnd[0], daughter1->PositionEnd[1], daughter1->PositionEnd[2]);
      TVector3 s2(daughter2->PositionStart[0], daughter2->PositionStart[1], daughter2->PositionStart[2]);
      TVector3 e2(daughter2->PositionEnd[0], daughter2->PositionEnd[1], daughter2->PositionEnd[2]);

      if (s1.X() < -900 || e1.X() < -900 || s2.X() < -900 || e2.X() < -900) {
        continue; // Skip if any endpoint invalid
      }

      std::vector<double> line1Params, line2Params;
      pdAnaUtils::ExtrapolateTrack(daughter1, line1Params, trackFitLength, true);
      pdAnaUtils::ExtrapolateTrack(daughter2, line2Params, trackFitLength, true);

      float distance;
      if (line1Params.size() >= 6 && line2Params.size() >= 6 &&
          line1Params[0] != -999.0 && line2Params[0] != -999.0) {
        TVector3 cp1, cp2;
        double closestApproach = pdAnaUtils::FindClosestPointsBetweenLines(line1Params, line2Params, cp1, cp2);
        distance = (closestApproach >= 0.0 && closestApproach < 1e6) ? static_cast<float>(closestApproach) : 1e6f;
      } else {
        // Fallback: min of four endpoint distances if track fit failed
        auto dist = [](const TVector3& a, const TVector3& b) { return (a - b).Mag(); };
        float d_ss = dist(s1, s2), d_se = dist(s1, e2), d_es = dist(e1, s2), d_ee = dist(e1, e2);
        distance = std::min({d_ss, d_se, d_es, d_ee});
      }

      if (distance > maxDaughterDistance) {
        continue; // Skip if tracks don't approach within threshold
      }

      // Create reconstructed vertex
      AnaAnnihilationVertexPD* reconstructedVertex = new AnaAnnihilationVertexPD();
      reconstructedVertex->OriginalDistance = distance;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->NParticles = 2;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);
      reconstructedVertex->ParentID = daughter1->ParentID;

      // Initialize vertex position and minimum distance to invalid values
      reconstructedVertex->Position[0] = -999.0;
      reconstructedVertex->Position[1] = -999.0;
      reconstructedVertex->Position[2] = -999.0;
      reconstructedVertex->MinimumDistance = -999.0;
      reconstructedVertex->Score = -999.0;

      // Call the position finder function to calculate vertex position and score
      positionFinder(reconstructedVertex, trackFitLength);

      // Calculate vertex direction (average of daughter directions)
      Float_t direction[3] = {daughter1->DirectionStart[0] + daughter2->DirectionStart[0],
                              daughter1->DirectionStart[1] + daughter2->DirectionStart[1],
                              daughter1->DirectionStart[2] + daughter2->DirectionStart[2]};

      Float_t norm = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
      direction[0] = direction[0] / norm;
      direction[1] = direction[1] / norm;
      direction[2] = direction[2] / norm;

      reconstructedVertex->Direction[0] = direction[0];
      reconstructedVertex->Direction[1] = direction[1];
      reconstructedVertex->Direction[2] = direction[2];

      // Compute the angle between the start directions of the two particles in degrees
      Float_t cosAngle = daughter1->DirectionStart[0]*daughter2->DirectionStart[0] +
                         daughter1->DirectionStart[1]*daughter2->DirectionStart[1] +
                         daughter1->DirectionStart[2]*daughter2->DirectionStart[2];
      // Clamp cosAngle to [-1, 1] to avoid numerical errors in acos
      if (cosAngle > 1.0) cosAngle = 1.0;
      if (cosAngle < -1.0) cosAngle = -1.0;
      reconstructedVertex->OpeningAngle = TMath::ACos(cosAngle) * 180.0 / TMath::Pi();

      // Calculate vertex momentum, energy, and angle with beam
      // Ensure particles have momentum calculated first
      reconstructedVertex->EnsureParticleMomentum();

      const double pionMass = 0.13957; // GeV
      TVector3 p1vec(daughter1->DirectionStart[0] * daughter1->Momentum,
                    daughter1->DirectionStart[1] * daughter1->Momentum,
                    daughter1->DirectionStart[2] * daughter1->Momentum);
      TVector3 p2vec(daughter2->DirectionStart[0] * daughter2->Momentum,
                    daughter2->DirectionStart[1] * daughter2->Momentum,
                    daughter2->DirectionStart[2] * daughter2->Momentum);
      TVector3 pTot = p1vec + p2vec;

      reconstructedVertex->Momentum[0] = pTot.X();
      reconstructedVertex->Momentum[1] = pTot.Y();
      reconstructedVertex->Momentum[2] = pTot.Z();

      double E1 = sqrt(daughter1->Momentum*daughter1->Momentum + pionMass*pionMass);
      double E2 = sqrt(daughter2->Momentum*daughter2->Momentum + pionMass*pionMass);
      reconstructedVertex->Energy = (E1 + E2) * 1000; // Convert to MeV

      // Calculate angle with beam direction
      TVector3 vtxDir(reconstructedVertex->Direction[0],
                     reconstructedVertex->Direction[1],
                     reconstructedVertex->Direction[2]);
      reconstructedVertex->AngleWithBeam = vtxDir.Angle(beamDir);

      // Create and fill true equivalent vertex using helper function
      reconstructedVertex->TrueEquivalentVertex = FillTrueEquivalentVertex(reconstructedVertex);

      // Set the generation and process
      reconstructedVertex->Generation = -999;
      reconstructedVertex->Process = -999;

      // Validate vertex before adding
      if (ValidateVertex(reconstructedVertex)) {
        reconstructedVertices.push_back(reconstructedVertex);
        // Mark this pair as used to prevent duplicate vertices
        usedPairs.insert(currentPair);
      } else {
        // Delete the vertex if validation failed
        delete reconstructedVertex;
      }
    }
  }

  // Filter vertices to ensure each particle belongs to at most one vertex
  std::vector<AnaAnnihilationVertexPD*> filteredVertices = FilterVerticesByScore(reconstructedVertices);

  return filteredVertices;
}

//***************************************************************
std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance) {
//***************************************************************

  // Check parameter to decide which algorithm to use
  int algorithmChoice = ND::params().GetParameterI("neutralKaonAnalysis.VertexFindingMethod");

  // Select the appropriate position finder function
  double (*positionFinder)(AnaVertexPD*, double) = nullptr;

  // 0 = Geometric, 1 = Fitted (TMinuit), 2 = Kalman
  switch(algorithmChoice) {
    case 0:
      positionFinder = FindVertexPositionGeometric;
      break;
    case 1:
      positionFinder = FindVertexPositionWithFit;
      break;
    case 2:
      positionFinder = FindVertexPositionKalman;
      break;
    default:
      std::cout << "WARNING: Unknown vertex finding method " << algorithmChoice
                << ", using Fitted method (1)" << std::endl;
      positionFinder = FindVertexPositionWithFit;
      break;
  }

  // Call common vertex creation function with selected position finder
  return CreateVerticesCommon(event, maxDaughterDistance, positionFinder);
}

//***************************************************************
void FindVertexPositionFit(AnaVertexPD* vertex){
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!part1 || !part2) {
    return;
  }

  // Get track fit length from parameters
  double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

  // Fit lines to both particles
  std::vector<double> line1Params, line2Params;
  pdAnaUtils::ExtrapolateTrack(part1, line1Params, trackFitLength, true); // Use start position for vertex particles
  pdAnaUtils::ExtrapolateTrack(part2, line2Params, trackFitLength, true); // Use start position for vertex particles

  // Check if both fits were successful
  if (line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    return;
  }

  // Find the closest points between the two fitted lines
  TVector3 closestPoint1, closestPoint2;
  double minDistance = pdAnaUtils::FindClosestPointsBetweenLines(line1Params, line2Params,
                                           closestPoint1, closestPoint2);

  // Set the vertex position to the midpoint between the closest points
  TVector3 vertexPosition = 0.5 * (closestPoint1 + closestPoint2);

  vertex->PositionFit[0] = vertexPosition.X();
  vertex->PositionFit[1] = vertexPosition.Y();
  vertex->PositionFit[2] = vertexPosition.Z();

  // Store closest points for downstream consumers (event display, etc.)
  vertex->ClosestPoint1Fit[0] = closestPoint1.X();
  vertex->ClosestPoint1Fit[1] = closestPoint1.Y();
  vertex->ClosestPoint1Fit[2] = closestPoint1.Z();
  vertex->ClosestPoint2Fit[0] = closestPoint2.X();
  vertex->ClosestPoint2Fit[1] = closestPoint2.Y();
  vertex->ClosestPoint2Fit[2] = closestPoint2.Z();

  // Store the minimum distance between the fitted lines
  vertex->MinimumDistanceFit = static_cast<Float_t>(minDistance);

  // Store the fitted line parameters for later use in event display
  vertex->FittedLineParams.clear();
  vertex->FittedLineParams.push_back(line1Params);
  vertex->FittedLineParams.push_back(line2Params);

}

//***************************************************************
void FindVertexPositionFit(AnaTrueEquivalentVertexPD* vertex){
//***************************************************************

  if (!vertex || vertex->TrueParticles.size() < 2) {
    return;
  }

  // Get the first two true particles from the vertex
  AnaTrueParticlePD* part1 = static_cast<AnaTrueParticlePD*>(vertex->TrueParticles[0]);
  AnaTrueParticlePD* part2 = static_cast<AnaTrueParticlePD*>(vertex->TrueParticles[1]);

  if (!part1 || !part2) {
    return;
  }

  // Get track fit length from parameters
  double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

  // Fit lines to both true particles
  std::vector<double> line1Params, line2Params;
  pdAnaUtils::ExtrapolateTrack(part1, line1Params, trackFitLength, true); // Use start position for vertex particles
  pdAnaUtils::ExtrapolateTrack(part2, line2Params, trackFitLength, true); // Use start position for vertex particles

  // Check if both fits were successful
  if (line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    return;
  }

  // Find the closest points between the two fitted lines
  TVector3 closestPoint1, closestPoint2;
  double minDistance = pdAnaUtils::FindClosestPointsBetweenLines(line1Params, line2Params,
                                                                 closestPoint1, closestPoint2);

  // Set the vertex position to the midpoint between the closest points
  TVector3 vertexPosition = 0.5 * (closestPoint1 + closestPoint2);

  vertex->PositionFit[0] = vertexPosition.X();
  vertex->PositionFit[1] = vertexPosition.Y();
  vertex->PositionFit[2] = vertexPosition.Z();

  vertex->MinimumDistance = static_cast<Float_t>(minDistance);

}

} // namespace pdAnnihilationUtils

