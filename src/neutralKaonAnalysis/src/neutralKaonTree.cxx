#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonAnalysisUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdDataClasses.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include "AnalysisUtils.hxx"
#include "ParticleId.hxx"
#include "HEPConstants.hxx"
#include "TH2F.h"
#include "TProfile.h"
#include "TH1F.h"
#include "CategoryManager.hxx"
#include "TMultiGraph.h"
#include "pdUtilsNeutralKaonSignalDiagnostics.hxx"
#include "pdMomReconstruction.hxx"
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>
#include <cstdio>

namespace {
  constexpr Float_t kHitTravelMinCm = 0.f;
  constexpr Float_t kHitTravelMaxCm = 500.f;

  const std::vector<Double_t>& GetHitTravelBinEdges() {
    static std::vector<Double_t> edges;
    if (!edges.empty()) {
      return edges;
    }

    // Fine granularity: 2 cm bins from 0 to 30 cm.
    for (Double_t x = kHitTravelMinCm; x < 30.; x += 2.) {
      edges.push_back(x);
    }

    // Coarser granularity: 5 cm bins from 30 to 500 cm.
    for (Double_t x = 30.; x <= kHitTravelMaxCm; x += 5.) {
      edges.push_back(x);
    }

    return edges;
  }

  TH2F* gK0DauHitDistVsTravel2D[2] = {nullptr, nullptr};
  TProfile* gK0DauHitDistVsTravelProf[2] = {nullptr, nullptr};
  TH2F* gK0SigDEdxVsRRByEndProc2D[4] = {nullptr, nullptr, nullptr, nullptr};
  const Int_t kTrackedEndProc[4] = {2, 9, 10, 11};

  Int_t EndProcToIndex(Int_t endProc) {
    for (Int_t i = 0; i < 4; ++i) {
      if (kTrackedEndProc[i] == endProc) return i;
    }
    return -1;
  }

  bool IsValidPos(const TVector3& p) {
    return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) &&
           p.X() > -900.f && p.Y() > -900.f && p.Z() > -900.f;
  }

  Float_t AngleBetweenDirections(TVector3 a, TVector3 b) {
    if (a.Mag2() <= 0.0 || b.Mag2() <= 0.0) return -999.0f;
    const double c = std::max(-1.0, std::min(1.0, a.Unit().Dot(b.Unit())));
    return static_cast<Float_t>(std::acos(c));
  }

  bool IsFinitePositive(const Float_t value) {
    return std::isfinite(value) && value > 0.f && value < 1e6f;
  }

  bool BuildMomentumVector(const Float_t magnitude, const TVector3& direction, TVector3& vector) {
    if (!IsFinitePositive(magnitude)) return false;
    if (!std::isfinite(direction.X()) || !std::isfinite(direction.Y()) || !std::isfinite(direction.Z())) return false;
    if (direction.Mag2() <= 0.) return false;
    vector = direction.Unit() * magnitude;
    return std::isfinite(vector.X()) && std::isfinite(vector.Y()) && std::isfinite(vector.Z()) && vector.Mag2() > 0.;
  }

  double SumCollectionPlaneVisibleEnergyMeV(const AnaParticlePD* part, int plane) {
    if (!part || plane < 0 || plane > 2) return 0.;
    double s = 0.;
    for (size_t i = 0; i < part->Hits[plane].size(); ++i) {
      const AnaHitPD& hit = part->Hits[plane][i];
      if (hit.dEdx <= 0.0f || hit.dEdx > 1000.0f || hit.dEdx == -999.0f) continue;
      double dx = 0.;
      if (hit.Pitch > 0.0f && hit.Pitch != -999.0f) {
        dx = hit.Pitch;
      } else {
        dx = 0.4792;
      }
      s += static_cast<double>(hit.dEdx) * dx;
    }
    return s;
  }

  double SumRecoDaughterSubtreesVisibleEnergyMeV(AnaParticlePD* vertexPion) {
    if (!vertexPion) return 0.;
    double total = 0.;
    for (size_t i = 0; i < vertexPion->Daughters.size(); ++i) {
      AnaParticlePD* root = static_cast<AnaParticlePD*>(vertexPion->Daughters[i]);
      if (!root) continue;
      std::vector<AnaParticlePD*> sub;
      pdMomReconstruction::CollectAllDescendants(root, sub);
      for (AnaParticlePD* node : sub) {
        total += SumCollectionPlaneVisibleEnergyMeV(node, 2);
      }
    }
    return total;
  }

  Float_t ComputeAlignmentAngle(const TVector3& creationPos, const TVector3& annihilationPos,
                                const TVector3& daughter1Momentum, const TVector3& daughter2Momentum) {
    if (!IsValidPos(creationPos) || !IsValidPos(annihilationPos)) return -999.f;

    const TVector3 neutralAxis = annihilationPos - creationPos;
    const TVector3 totalMomentum = daughter1Momentum + daughter2Momentum;

    if (neutralAxis.Mag2() <= 0. || totalMomentum.Mag2() <= 0.) return -999.f;

    const double c = neutralAxis.Dot(totalMomentum) / (neutralAxis.Mag() * totalMomentum.Mag());
    if (!std::isfinite(c)) return -999.f;
    const double cClamped = std::max(-1.0, std::min(1.0, c));
    return static_cast<Float_t>(std::acos(cClamped));
  }

  void EnsureHitDistanceProfiles(OutputManager& output) {
    (void)output;
    if (gK0DauHitDistVsTravel2D[0] && gK0DauHitDistVsTravel2D[1] &&
        gK0DauHitDistVsTravelProf[0] && gK0DauHitDistVsTravelProf[1]) {
      return;
    }
    for (Int_t i = 0; i < 2; ++i) {
      const Int_t dauIdx = i + 1;
      const std::vector<Double_t>& travelEdges = GetHitTravelBinEdges();
      const Int_t nTravelBins = static_cast<Int_t>(travelEdges.size()) - 1;
      if (!gK0DauHitDistVsTravel2D[i]) {
        gK0DauHitDistVsTravel2D[i] = new TH2F(
            Form("h_k0dau%d_hitDistToTrueLine_vs_travel_2d", dauIdx),
            Form("K0 daughter %d: hit distance to true line vs traveled distance;Traveled distance from reco start [cm];Hit-to-true-line distance [cm]", dauIdx),
            nTravelBins, travelEdges.data(),
            200, 0.f, 100.f);
        gK0DauHitDistVsTravel2D[i]->SetDirectory(nullptr);
      }

      if (!gK0DauHitDistVsTravelProf[i]) {
        gK0DauHitDistVsTravelProf[i] = new TProfile(
            Form("p_k0dau%d_hitDistToTrueLine_vs_travel", dauIdx),
            Form("K0 daughter %d: mean hit distance to true line vs traveled distance;Traveled distance from reco start [cm];Mean hit-to-true-line distance [cm]", dauIdx),
            nTravelBins, travelEdges.data());
        gK0DauHitDistVsTravelProf[i]->SetDirectory(nullptr);
      }
    }
  }

  void EnsureDEdxResidualRangeByEndProcHistograms(OutputManager& output) {
    (void)output;
    for (Int_t i = 0; i < 4; ++i) {
      if (!gK0SigDEdxVsRRByEndProc2D[i]) {
        gK0SigDEdxVsRRByEndProc2D[i] = new TH2F(
            Form("h_k0sig_dEdx_vs_rr_endproc%d_2d", kTrackedEndProc[i]),
            Form("Signal K0 daughters: dE/dx vs residual range (true end proc %d);Residual range [cm];dE/dx [MeV/cm]",
                 kTrackedEndProc[i]),
            150, 0.f, 150.f,
            300, 0.f, 30.f);
        gK0SigDEdxVsRRByEndProc2D[i]->SetDirectory(nullptr);
      }
    }
  }

  void FillHitDistanceToTrueLineProfiles(OutputManager& output, const AnaParticlePD* recoParticle,
                                         const AnaTrueParticlePD* trueParticle, Int_t daughterIndex) {
    if (!recoParticle || !trueParticle || daughterIndex < 0 || daughterIndex > 1) {
      return;
    }

    EnsureHitDistanceProfiles(output);
    if (!gK0DauHitDistVsTravel2D[daughterIndex] || !gK0DauHitDistVsTravelProf[daughterIndex]) {
      return;
    }

    const TVector3 recoStart(recoParticle->PositionStart[0], recoParticle->PositionStart[1], recoParticle->PositionStart[2]);
    const TVector3 trueStart(trueParticle->Position[0], trueParticle->Position[1], trueParticle->Position[2]);
    TVector3 trueDir(trueParticle->Direction[0], trueParticle->Direction[1], trueParticle->Direction[2]);

    if (!IsValidPos(recoStart) || !IsValidPos(trueStart)) {
      return;
    }
    if (trueDir.Mag2() <= 0.f) {
      return;
    }
    trueDir = trueDir.Unit();

    const std::vector<AnaHitPD>& hits = recoParticle->Hits[2];
    for (size_t ihit = 0; ihit < hits.size(); ++ihit) {
      const TVector3 hitPos = hits[ihit].Position;
      if (!IsValidPos(hitPos)) {
        continue;
      }

      const Float_t traveledDistance = static_cast<Float_t>((hitPos - recoStart).Mag());
      if (traveledDistance < kHitTravelMinCm || traveledDistance >= kHitTravelMaxCm) {
        continue;
      }

      const Float_t hitToTrueLineDistance = static_cast<Float_t>(((hitPos - trueStart).Cross(trueDir)).Mag());
      gK0DauHitDistVsTravel2D[daughterIndex]->Fill(traveledDistance, hitToTrueLineDistance);
      gK0DauHitDistVsTravelProf[daughterIndex]->Fill(traveledDistance, hitToTrueLineDistance);
    }
  }

  void FillDEdxVsResidualRangeByEndProc(OutputManager& output, const AnaParticlePD* recoParticle,
                                        const AnaTrueParticlePD* trueParticle) {
    if (!recoParticle || !trueParticle) return;
    const Int_t procIdx = EndProcToIndex(trueParticle->ProcessEnd);
    if (procIdx < 0 || procIdx > 3) return;

    EnsureDEdxResidualRangeByEndProcHistograms(output);
    if (!gK0SigDEdxVsRRByEndProc2D[procIdx]) return;

    const std::vector<AnaHitPD>& hits = recoParticle->Hits[2];
    for (size_t ihit = 0; ihit < hits.size(); ++ihit) {
      const AnaHitPD& hit = hits[ihit];
      const Float_t rr = hit.ResidualRange;
      const Float_t dedx = (hit.dEdx_calib > 0.0f && hit.dEdx_calib != -999.0f && std::isfinite(hit.dEdx_calib))
                               ? hit.dEdx_calib
                               : hit.dEdx;
      if (!(rr >= 0.0f) || rr == -999.0f || !std::isfinite(rr)) continue;
      if (!(dedx > 0.0f) || dedx == -999.0f || dedx >= 1000.0f || !std::isfinite(dedx)) continue;

      gK0SigDEdxVsRRByEndProc2D[procIdx]->Fill(rr, dedx);
    }
  }

  Int_t CountRecoObjectsForTrueDaughters(const AnaTrueParticlePD* trueParticle, const AnaEventB& event) {
    if (!trueParticle) {
      return -1;
    }

    Int_t nRecoDaughters = 0;
    for (size_t idau = 0; idau < trueParticle->Daughters.size(); ++idau) {
      const Int_t trueDaughterID = trueParticle->Daughters[idau];

      bool hasRecoMatch = false;
      for (Int_t ip = 0; ip < event.nParticles; ++ip) {
        const AnaParticlePD* recoParticle = static_cast<const AnaParticlePD*>(event.Particles[ip]);
        if (!recoParticle || !recoParticle->TrueObject) {
          continue;
        }
        if (recoParticle->TrueObject->ID == trueDaughterID) {
          hasRecoMatch = true;
          break;
        }
      }

      if (hasRecoMatch) {
        ++nRecoDaughters;
      }
    }
    return nRecoDaughters;
  }

}

//********************************************************************
void neutralKaonTree::WriteHitDistanceProfiles(OutputManager& output){
//********************************************************************
  // EnsureDEdxResidualRangeByEndProcHistograms(output);
  // for (Int_t i = 0; i < 2; ++i) {
  //   if (gK0DauHitDistVsTravel2D[i]) {
  //     gK0DauHitDistVsTravel2D[i]->Write("", TObject::kOverwrite);
  //   }
  //   if (gK0DauHitDistVsTravelProf[i]) {
  //     gK0DauHitDistVsTravelProf[i]->Write("", TObject::kOverwrite);
  //   }
  // }
  // for (Int_t i = 0; i < 4; ++i) {
  //   if (gK0SigDEdxVsRRByEndProc2D[i]) {
  //     gK0SigDEdxVsRRByEndProc2D[i]->Write("", TObject::kOverwrite);
  //   }
  // }
  neutralKaonTreeDiagnostics::WriteSignalPionDedxDiagnostics(output);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Particle(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVF(output, k0lengthpandora, "Neutral length using annihilation Pandora position", nk0, nmax);
  AddVarMaxSizeVF(output, k0lengthfit, "Neutral length using annihilation Fit position", nk0, nmax);
  AddVarMaxSizeVF(output, k0truelength, "True K0 length using true creation and decay vertices", nk0, nmax);
  AddVarMaxSizeVF(output, k0alignmentpandora,
                   "K0 alignment (Pandora): angle (rad) between creation→annihilation(Pandora) and vertex Σp (Pandora dirs)", nk0,
                   nmax);
  AddVarMaxSizeVF(output, k0alignmentfit,
                   "K0 alignment (Fit): angle (rad) between creation→annihilation(Fit) and vertex Σp (fit dirs)", nk0,
                   nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0AlignVariants(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVF(output, k0al_alltrue, "Alignment with all true quantities", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_allreco, "Alignment with all reco quantities (fit reco geometry/directions)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_cvreco, "Alignment with reco creation only (all else true)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_avreco, "Alignment with reco annihilation only (all else true)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_vtxreco, "Alignment with reco creation+annihilation only (all else true)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d1magreco, "Alignment with daughter1 reco momentum magnitude only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d2magreco, "Alignment with daughter2 reco momentum magnitude only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d1dirreco, "Alignment with daughter1 reco momentum direction only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d2dirreco, "Alignment with daughter2 reco momentum direction only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d1preco, "Alignment with daughter1 reco momentum (direction+magnitude)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d2preco, "Alignment with daughter2 reco momentum (direction+magnitude)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d12magreco, "Alignment with both daughters reco momentum magnitudes only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d12dirreco, "Alignment with both daughters reco momentum directions only", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_d12preco, "Alignment with both daughters reco momentum (direction+magnitude)", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_allreco_cvtrue, "Alignment with all reco except creation vertex true", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_allreco_avtrue, "Alignment with all reco except annihilation vertex true", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_allreco_d1true, "Alignment with all reco except daughter1 momentum true", nk0, nmax);
  AddVarMaxSizeVF(output, k0al_allreco_d2true, "Alignment with all reco except daughter2 momentum true", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Parent(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVI(output, k0partruepdg, "True PDG of K0 parent", nk0, nmax);
  AddVarMaxSizeVF(output, k0partrueendmom, "True end momentum of K0 parent", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruelength, "True track length of K0 parent", nk0, nmax);
  AddVarMaxSizeVF(output, k0parrecolength, "Reco track length of K0 parent", nk0, nmax);
  AddVarMaxSizeVI(output, k0parndau, "Reco number of daughters of K0 parent", nk0, nmax);
  AddVarMaxSizeVI(output, k0parisbeam, "1 if reco parent IsPandora (beam), 0 otherwise", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0CreationVtx(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVF(output, k0cvtxpandoraresidual,
                  "Creation-vertex residual |Pandora(raw parent end) - true creation|", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxfitresidual,
                  "Creation-vertex residual |Fit(projected parent end) - true creation|", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxpandoraresidualx, "Creation Pandora residual x_reco - x_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxpandoraresidualy, "Creation Pandora residual y_reco - y_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxpandoraresidualz, "Creation Pandora residual z_reco - z_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxfitresidualx, "Creation Fit residual x_reco - x_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxfitresidualy, "Creation Fit residual y_reco - y_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0cvtxfitresidualz, "Creation Fit residual z_reco - z_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0protonmomentumreco, "Reco momentum magnitude of creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVF(output, k0protonmomentumtrue, "True momentum magnitude of creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVI(output, k0protontruepdg, "True PDG of creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVF(output, k0protonchi2ndfproton, "Chi2PID(2212)/npts creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVF(output, k0protontruelength, "True track length of creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVF(output, k0protonrecolength, "Reco track length of creation-vertex second particle", nk0, nmax);
  AddVarMaxSizeVI(output, k0hasproton, "1 if creation-vertex second particle is assigned and valid, 0 otherwise", nk0, nmax);
  AddVarMaxSizeVI(output, k0creationdegeneracy, "K0 creation-vertex degeneracy (Reco)", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax){

  AddVarMaxSizeVI(output, k0nvtxbeforefiltering, "Number of annihilation vertices before overlap filtering", nk0, nmax);
  AddVarMaxSizeVI(output, k0nvtxafterfiltering, "Number of annihilation vertices after overlap filtering", nk0, nmax);

  //Vertex system variables
  AddVarMaxSize3MF(output, k0vtxtruepos, "K0 vertex true position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxoriginaldistance, "K0 vertex original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueoriginaldistance, "K0 vertex true original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdegeneracy, "K0 annihilation-vertex degeneracy (Reco)", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpandorapos, "K0 vertex Pandora position", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxfitpos, "K0 vertex fitted position (geometric/TMinuit/Kalman)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandorax, "K0 vertex Pandora x position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoray, "K0 vertex Pandora y position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoraz, "K0 vertex Pandora z position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitx, "K0 vertex Fit x position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfity, "K0 vertex Fit y position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitz, "K0 vertex Fit z position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoraresidual, "Vertex true-reco distance (Pandora position)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitresidual, "Vertex true-reco distance (Fit position)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoraresidualx, "Vertex Pandora residual x_reco - x_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoraresidualy, "Vertex Pandora residual y_reco - y_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpandoraresidualz, "Vertex Pandora residual z_reco - z_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitresidualx, "Vertex Fit residual x_reco - x_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitresidualy, "Vertex Fit residual y_reco - y_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxfitresidualz, "Vertex Fit residual z_reco - z_true", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxmomentum, "K0 annihilation-vertex total momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxinvariantmass, "K0 annihilation-vertex invariant mass (pion hypothesis)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxmomentumpandora, "K0 vertex momentum using Pandora directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxinvariantmasspandora, "K0 vertex invariant mass using Pandora directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxopeninganglepandora, "K0 vertex opening angle using Pandora daughter directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxmomentumfit, "K0 vertex momentum using fit directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxinvariantmassfit, "K0 vertex invariant mass using fit directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxopeninganglefit, "K0 vertex opening angle using fit daughter directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxresultantmomentumreco, "Reco resultant momentum magnitude of annihilation vertex", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxresultantmomentumtrue, "True resultant momentum magnitude from true daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxjointk0smomused, "1 if joint K0s momentum was applied for this vertex", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointsigmap1, "Joint-fit sigma(p1) from TLE curvature at marginal best [GeV/c]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointsigmap2, "Joint-fit sigma(p2) from TLE curvature at marginal best [GeV/c]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointsigmam, "Joint-fit event sigma_m used in mass penalty [GeV]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointdmdp1, "Joint-fit dm/dp1 at marginal (p1_0,p2_0)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointdmdp2, "Joint-fit dm/dp2 at marginal (p1_0,p2_0)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointrconstraint, "Joint post-fit kinematic constraint dominance ratio R", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxjointdeltachidedx, "Joint post-fit dE/dx template chi2 degradation (or chi2 at joint)", nk0,
                  nmax);
  AddVarMaxSizeVI(output, k0vtxjointdebugclass, "Joint post-fit debug classification (see enum in neutralKaonTree.hxx)",
                  nk0, nmax);
  AddNeutralKaonVariables_K0VtxDaughters(output, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0VtxDaughters(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVF(output, k0vtxdau1momentumreco, "Reco momentum magnitude of daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2momentumreco, "Reco momentum magnitude of daughter 2", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1mommethod, "Momentum assignment method enum for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2mommethod, "Momentum assignment method enum for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1extchi2ndf, "Momentum-assignment free-range fit log-likelihood for daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2extchi2ndf, "Momentum-assignment free-range fit log-likelihood for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1dedxdrift, "Mean dEdx bias (actual - expected from fit PDF mode) for daughter 1 (MeV/cm)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2dedxdrift, "Mean dEdx bias (actual - expected from fit PDF mode) for daughter 2 (MeV/cm)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1dedxsigma, "Sigma of dEdx bias Gaussian fit for daughter 1 (MeV/cm)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2dedxsigma, "Sigma of dEdx bias Gaussian fit for daughter 2 (MeV/cm)", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1dedxfitok, "Gaussian dEdx bias fit success flag for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2dedxfitok, "Gaussian dEdx bias fit success flag for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1momentumtrue, "True momentum magnitude of daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2momentumtrue, "True momentum magnitude of daughter 2", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1trueendproc, "True end process enum for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2trueendproc, "True end process enum for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1trueendmom, "True end momentum for daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2trueendmom, "True end momentum for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truestartmom, "True start momentum for daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truestartmom, "True start momentum for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truestartx, "True start x for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truestarty, "True start y for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truestartz, "True start z for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truestartx, "True start x for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truestarty, "True start y for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truestartz, "True start z for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1truendau, "True number of daughters for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2truendau, "True number of daughters for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truelength, "True track length for daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truelength, "True track length for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1recolength, "Reco track length for daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2recolength, "Reco track length for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1recostartx, "Reco start x for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1recostarty, "Reco start y for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1recostartz, "Reco start z for daughter 1 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2recostartx, "Reco start x for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2recostarty, "Reco start y for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2recostartz, "Reco start z for daughter 2 [cm]", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1nhitsreco, "Reco collection-plane hit count for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2nhitsreco, "Reco collection-plane hit count for daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1protonchi2ndf, "Chi2PID(2212)/npts vertex daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2protonchi2ndf, "Chi2PID(2212)/npts vertex daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1pionchi2ndf, "Chi2PID(211)/npts vertex daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2pionchi2ndf, "Chi2PID(211)/npts vertex daughter 2", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1kaonchi2ndf, "Chi2PID(321)/npts vertex daughter 1", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2kaonchi2ndf, "Chi2PID(321)/npts vertex daughter 2", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1truepdg, "True PDG for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2truepdg, "True PDG for daughter 2", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1ndaughtersreco, "Reco number of daughters for daughter 1", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2ndaughtersreco, "Reco number of daughters for daughter 2", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau1nrecodau, "Number of true daughters of daughter 1 with reco match", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdau2nrecodau, "Number of true daughters of daughter 2 with reco match", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1recovisiblee,
                  "Visible calo energy [GeV] in reco daughter subtree of vertex pion 1 (coll. plane dE/dx*dx)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2recovisiblee,
                  "Visible calo energy [GeV] in reco daughter subtree of vertex pion 2 (coll. plane dE/dx*dx)", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdaughtersrecovisiblee,
                  "Visible calo energy [GeV] in reco daughter subtrees of both vertex pions (coll. plane dE/dx*dx)", nk0,
                  nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Kinematics(OutputManager& output, UInt_t nmax){
//********************************************************************
  // True K0 kinematics
  AddVarMaxSizeVF(output, k0truecreationx, "True K0 creation position x [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truecreationy, "True K0 creation position y [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truecreationz, "True K0 creation position z [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0trueannihilationx, "True K0 annihilation position x [cm]", nk0, nmax);
    AddVarMaxSizeVF(output, k0trueannihilationy, "True K0 annihilation position y [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0trueannihilationz, "True K0 annihilation position z [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truemomentum, "True K0 momentum magnitude [GeV/c]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truedirectionx, "True K0 direction (normalized) x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truedirectiony, "True K0 direction (normalized) y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0truedirectionz, "True K0 direction (normalized) z [dimensionless]", nk0, nmax);

  // Reco K0 creation position: Pandora variant
  AddVarMaxSizeVF(output, k0creationpandorax, "Reco K0 creation Pandora x [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0creationpandoray, "Reco K0 creation Pandora y [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0creationpandoraz, "Reco K0 creation Pandora z [cm]", nk0, nmax);

  // Reco K0 creation position: Fit variant
  AddVarMaxSizeVF(output, k0creationfitx, "Reco K0 creation Fit x [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0creationfity, "Reco K0 creation Fit y [cm]", nk0, nmax);
  AddVarMaxSizeVF(output, k0creationfitz, "Reco K0 creation Fit z [cm]", nk0, nmax);

  // Reco K0 direction: Pandora variant
  AddVarMaxSizeVF(output, k0directionpandorax, "Reco K0 direction Pandora x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0directionpandoray, "Reco K0 direction Pandora y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0directionpandoraz, "Reco K0 direction Pandora z [dimensionless]", nk0, nmax);

  // Reco K0 direction: Fit variant
  AddVarMaxSizeVF(output, k0directionfitx, "Reco K0 direction Fit x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0directionfity, "Reco K0 direction Fit y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0directionfitz, "Reco K0 direction Fit z [dimensionless]", nk0, nmax);

  // True daughter directions
  AddVarMaxSizeVF(output, k0vtxdau1truedirectionx, "True daughter 1 direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truedirectiony, "True daughter 1 direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1truedirectionz, "True daughter 1 direction z [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truedirectionx, "True daughter 2 direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truedirectiony, "True daughter 2 direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2truedirectionz, "True daughter 2 direction z [dimensionless]", nk0, nmax);

  // Reco daughter directions: Pandora variant
  AddVarMaxSizeVF(output, k0vtxdau1directionpandorax, "Reco daughter 1 Pandora direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1directionpandoray, "Reco daughter 1 Pandora direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1directionpandoraz, "Reco daughter 1 Pandora direction z [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionpandorax, "Reco daughter 2 Pandora direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionpandoray, "Reco daughter 2 Pandora direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionpandoraz, "Reco daughter 2 Pandora direction z [dimensionless]", nk0, nmax);

  // Reco daughter directions: Fit variant
  AddVarMaxSizeVF(output, k0vtxdau1directionfitx, "Reco daughter 1 Fit direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1directionfity, "Reco daughter 1 Fit direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau1directionfitz, "Reco daughter 1 Fit direction z [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionfitx, "Reco daughter 2 Fit direction x [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionfity, "Reco daughter 2 Fit direction y [dimensionless]", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxdau2directionfitz, "Reco daughter 2 Fit direction z [dimensionless]", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                               Int_t nVerticesBeforeFiltering, Int_t nVerticesAfterFiltering,
                                               AnaBeamB* beam, size_t neutralCandidateIndex){
  (void)beam;
  output.FillVectorVar(k0nvtxbeforefiltering, nVerticesBeforeFiltering);
  output.FillVectorVar(k0nvtxafterfiltering, nVerticesAfterFiltering);
  neutralKaonTree::FillNeutralKaonVariables_K0Particle(output, candidate);
  neutralKaonTree::FillNeutralKaonVariables_K0AlignVariants(output, candidate);
  neutralKaonTree::FillNeutralKaonVariables_K0Parent(output, candidate);
  neutralKaonTree::FillNeutralKaonVariables_K0CreationVtx(output, candidate);
  neutralKaonTree::FillNeutralKaonVariables_K0Kinematics(output, candidate);

  if(candidate){
    AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
    const Int_t excludedParentUniqueID = candidate->Parent ? candidate->Parent->UniqueID : -1;
    neutralKaonTree::FillNeutralKaonVariables_K0vtx(output, vertex, event, excludedParentUniqueID,
                                                    candidate->CreationVertex);

    // neutralKaonTreeDiagnostics::MaybeAccumulateSignalPionDedxMultiGraphs(candidate, event, neutralCandidateIndex);

    AnaTrueParticlePD* sigTrueK0 = neutralKaonAnaUtils::GetSignalTrueParent(candidate, event);
    if (vertex && sigTrueK0) {
      neutralKaonTreeDiagnostics::ResetSignalK0DiagnosticsIfNewEvent(event);
      if (neutralKaonTreeDiagnostics::RegisterSignalTrueK0Id(sigTrueK0->ID)) {
        if (vertex->Particles.size() > 0 && vertex->Particles[0]) {
          AnaParticlePD* recoDau1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
          AnaTrueParticlePD* trueDau1 = recoDau1 ? static_cast<AnaTrueParticlePD*>(recoDau1->TrueObject) : nullptr;
          // FillHitDistanceToTrueLineProfiles(output, recoDau1, trueDau1, 0);
          // FillDEdxVsResidualRangeByEndProc(output, recoDau1, trueDau1);
        }
        if (vertex->Particles.size() > 1 && vertex->Particles[1]) {
          AnaParticlePD* recoDau2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
          AnaTrueParticlePD* trueDau2 = recoDau2 ? static_cast<AnaTrueParticlePD*>(recoDau2->TrueObject) : nullptr;
          // FillHitDistanceToTrueLineProfiles(output, recoDau2, trueDau2, 1);
          // FillDEdxVsResidualRangeByEndProc(output, recoDau2, trueDau2);
        }
      }
    }
  }
}
//********************************************************************

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Particle(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
  Float_t lengthPandora = -999.0f;
  Float_t lengthFit = -999.0f;
  Float_t trueLength = -999.0f;
  Float_t alignmentPandora = -999.0f;
  Float_t alignmentFit = -999.0f;

  if (candidate) {
    lengthPandora = candidate->LengthPandora;
    lengthFit = candidate->LengthFit;
    AnaTrueParticlePD* trueK0 = candidate->TrueObject ? static_cast<AnaTrueParticlePD*>(candidate->TrueObject) : nullptr;
    if (trueK0) {
      const TVector3 trueStart(trueK0->Position[0], trueK0->Position[1], trueK0->Position[2]);
      const TVector3 trueEnd(trueK0->PositionEnd[0], trueK0->PositionEnd[1], trueK0->PositionEnd[2]);
      if (IsValidPos(trueStart) && IsValidPos(trueEnd)) {
        trueLength = static_cast<Float_t>((trueEnd - trueStart).Mag());
      }
    }
    alignmentPandora = candidate->AlignmentPandora;
    alignmentFit = candidate->AlignmentFit;
  }

  output.FillVectorVar(k0lengthpandora, lengthPandora);
  output.FillVectorVar(k0lengthfit, lengthFit);
  output.FillVectorVar(k0truelength, trueLength);
  output.FillVectorVar(k0alignmentpandora, alignmentPandora);
  output.FillVectorVar(k0alignmentfit, alignmentFit);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0AlignVariants(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
  Float_t k0al_alltrue_val = -999.f;
  Float_t k0al_allreco_val = -999.f;
  Float_t k0al_cvreco_val = -999.f;
  Float_t k0al_avreco_val = -999.f;
  Float_t k0al_vtxreco_val = -999.f;
  Float_t k0al_d1magreco_val = -999.f;
  Float_t k0al_d2magreco_val = -999.f;
  Float_t k0al_d1dirreco_val = -999.f;
  Float_t k0al_d2dirreco_val = -999.f;
  Float_t k0al_d1preco_val = -999.f;
  Float_t k0al_d2preco_val = -999.f;
  Float_t k0al_d12magreco_val = -999.f;
  Float_t k0al_d12dirreco_val = -999.f;
  Float_t k0al_d12preco_val = -999.f;
  Float_t k0al_allreco_cvtrue_val = -999.f;
  Float_t k0al_allreco_avtrue_val = -999.f;
  Float_t k0al_allreco_d1true_val = -999.f;
  Float_t k0al_allreco_d2true_val = -999.f;

  if (candidate) {
    AnaCreationVertexPD* creationVtx = candidate->CreationVertex;
    AnaParticlePD* recoParent = nullptr;
    if (creationVtx && creationVtx->BeamParticle) {
      recoParent = creationVtx->BeamParticle;
    } else if (candidate->Parent) {
      recoParent = static_cast<AnaParticlePD*>(candidate->Parent);
    }
    AnaAnnihilationVertexPD* annihilationVtx = candidate->AnnihilationVertex;

    AnaTrueParticlePD* trueK0 = candidate->TrueObject ? static_cast<AnaTrueParticlePD*>(candidate->TrueObject) : nullptr;
    AnaTrueParticlePD* trueParent =
        (recoParent && recoParent->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoParent->TrueObject) : nullptr;

    TVector3 trueCreation(-999., -999., -999.);
    bool hasTrueCreation = false;
    if (trueK0) {
      trueCreation.SetXYZ(trueK0->Position[0], trueK0->Position[1], trueK0->Position[2]);
      hasTrueCreation = IsValidPos(trueCreation);
    }
    if (!hasTrueCreation && trueParent) {
      trueCreation.SetXYZ(trueParent->PositionEnd[0], trueParent->PositionEnd[1], trueParent->PositionEnd[2]);
      hasTrueCreation = IsValidPos(trueCreation);
    }

    TVector3 trueAnnihilation(-999., -999., -999.);
    bool hasTrueAnnihilation = false;
    if (trueK0) {
      trueAnnihilation.SetXYZ(trueK0->PositionEnd[0], trueK0->PositionEnd[1], trueK0->PositionEnd[2]);
      hasTrueAnnihilation = IsValidPos(trueAnnihilation);
    }

    TVector3 recoCreation(-999., -999., -999.);
    bool hasRecoCreation = false;
    if (creationVtx) {
      recoCreation.SetXYZ(creationVtx->Position[0], creationVtx->Position[1], creationVtx->Position[2]);
      hasRecoCreation = IsValidPos(recoCreation);
    }

    TVector3 recoAnnihilation(-999., -999., -999.);
    bool hasRecoAnnihilation = false;
    if (annihilationVtx) {
      recoAnnihilation.SetXYZ(annihilationVtx->PositionFit[0], annihilationVtx->PositionFit[1], annihilationVtx->PositionFit[2]);
      hasRecoAnnihilation = IsValidPos(recoAnnihilation);
    }

    AnaParticlePD* recoDau1 = nullptr;
    AnaParticlePD* recoDau2 = nullptr;
    AnaTrueParticlePD* trueDau1 = nullptr;
    AnaTrueParticlePD* trueDau2 = nullptr;
    if (annihilationVtx && annihilationVtx->Particles.size() >= 2) {
      recoDau1 = static_cast<AnaParticlePD*>(annihilationVtx->Particles[0]);
      recoDau2 = static_cast<AnaParticlePD*>(annihilationVtx->Particles[1]);
      trueDau1 = (recoDau1 && recoDau1->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoDau1->TrueObject) : nullptr;
      trueDau2 = (recoDau2 && recoDau2->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoDau2->TrueObject) : nullptr;
    }

    Float_t d1TrueMag = -999.f, d2TrueMag = -999.f;
    Float_t d1RecoMag = -999.f, d2RecoMag = -999.f;
    TVector3 d1TrueDir(0., 0., 0.), d2TrueDir(0., 0., 0.);
    TVector3 d1RecoDir(0., 0., 0.), d2RecoDir(0., 0., 0.);
    bool hasD1TrueMag = false, hasD2TrueMag = false;
    bool hasD1RecoMag = false, hasD2RecoMag = false;
    bool hasD1TrueDir = false, hasD2TrueDir = false;
    bool hasD1RecoDir = false, hasD2RecoDir = false;

    if (trueDau1) {
      d1TrueMag = trueDau1->Momentum;
      hasD1TrueMag = IsFinitePositive(d1TrueMag);
      d1TrueDir.SetXYZ(trueDau1->Direction[0], trueDau1->Direction[1], trueDau1->Direction[2]);
      hasD1TrueDir = (d1TrueDir.Mag2() > 0. && std::isfinite(d1TrueDir.X()) && std::isfinite(d1TrueDir.Y()) && std::isfinite(d1TrueDir.Z()));
    }
    if (trueDau2) {
      d2TrueMag = trueDau2->Momentum;
      hasD2TrueMag = IsFinitePositive(d2TrueMag);
      d2TrueDir.SetXYZ(trueDau2->Direction[0], trueDau2->Direction[1], trueDau2->Direction[2]);
      hasD2TrueDir = (d2TrueDir.Mag2() > 0. && std::isfinite(d2TrueDir.X()) && std::isfinite(d2TrueDir.Y()) && std::isfinite(d2TrueDir.Z()));
    }
    const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
    const double trackFitDistanceFromStart = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");

    if (recoDau1) {
      d1RecoMag = recoDau1->Momentum;
      hasD1RecoMag = IsFinitePositive(d1RecoMag);
      d1RecoDir.SetXYZ(recoDau1->DirectionStart[0], recoDau1->DirectionStart[1], recoDau1->DirectionStart[2]);
      hasD1RecoDir = (d1RecoDir.Mag2() > 0. && std::isfinite(d1RecoDir.X()) && std::isfinite(d1RecoDir.Y()) && std::isfinite(d1RecoDir.Z()));

      std::vector<double> fit1;
      pdAnaUtils::ExtrapolateTrack(recoDau1, fit1, trackFitLength, true, trackFitDistanceFromStart);
      const bool fit1Valid =
          (fit1.size() >= 6 && fit1[3] > -900.0 && fit1[4] > -900.0 && fit1[5] > -900.0 && std::isfinite(fit1[3]) &&
           std::isfinite(fit1[4]) && std::isfinite(fit1[5]));
      if (fit1Valid) {
        d1RecoDir.SetXYZ(fit1[3], fit1[4], fit1[5]);
        hasD1RecoDir = true;
      }
    }
    if (recoDau2) {
      d2RecoMag = recoDau2->Momentum;
      hasD2RecoMag = IsFinitePositive(d2RecoMag);
      d2RecoDir.SetXYZ(recoDau2->DirectionStart[0], recoDau2->DirectionStart[1], recoDau2->DirectionStart[2]);
      hasD2RecoDir = (d2RecoDir.Mag2() > 0. && std::isfinite(d2RecoDir.X()) && std::isfinite(d2RecoDir.Y()) && std::isfinite(d2RecoDir.Z()));

      std::vector<double> fit2;
      pdAnaUtils::ExtrapolateTrack(recoDau2, fit2, trackFitLength, true, trackFitDistanceFromStart);
      const bool fit2Valid =
          (fit2.size() >= 6 && fit2[3] > -900.0 && fit2[4] > -900.0 && fit2[5] > -900.0 && std::isfinite(fit2[3]) &&
           std::isfinite(fit2[4]) && std::isfinite(fit2[5]));
      if (fit2Valid) {
        d2RecoDir.SetXYZ(fit2[3], fit2[4], fit2[5]);
        hasD2RecoDir = true;
      }
    }

    auto buildDaughterMomentum = [&](bool useRecoMag, bool useRecoDir,
                                     Float_t trueMag, bool hasTrueMag, const TVector3& trueDir, bool hasTrueDir,
                                     Float_t recoMag, bool hasRecoMag, const TVector3& recoDir, bool hasRecoDir,
                                     TVector3& out) {
      const Float_t mag = useRecoMag ? recoMag : trueMag;
      const TVector3 dir = useRecoDir ? recoDir : trueDir;
      const bool magOk = useRecoMag ? hasRecoMag : hasTrueMag;
      const bool dirOk = useRecoDir ? hasRecoDir : hasTrueDir;
      if (!magOk || !dirOk) return false;
      return BuildMomentumVector(mag, dir, out);
    };

    auto computeVariant = [&](Float_t& out,
                  bool useRecoCreation, bool useRecoAnnihilation,
                  bool useD1RecoMag, bool useD1RecoDir,
                  bool useD2RecoMag, bool useD2RecoDir) {
      const bool cOk = useRecoCreation ? hasRecoCreation : hasTrueCreation;
      const bool aOk = useRecoAnnihilation ? hasRecoAnnihilation : hasTrueAnnihilation;
      if (!cOk || !aOk) return;

      const TVector3 cPos = useRecoCreation ? recoCreation : trueCreation;
      const TVector3 aPos = useRecoAnnihilation ? recoAnnihilation : trueAnnihilation;

      TVector3 p1(0., 0., 0.), p2(0., 0., 0.);
      if (!buildDaughterMomentum(useD1RecoMag, useD1RecoDir,
                                 d1TrueMag, hasD1TrueMag, d1TrueDir, hasD1TrueDir,
                                 d1RecoMag, hasD1RecoMag, d1RecoDir, hasD1RecoDir,
                                 p1)) return;
      if (!buildDaughterMomentum(useD2RecoMag, useD2RecoDir,
                                 d2TrueMag, hasD2TrueMag, d2TrueDir, hasD2TrueDir,
                                 d2RecoMag, hasD2RecoMag, d2RecoDir, hasD2RecoDir,
                                 p2)) return;

      out = ComputeAlignmentAngle(cPos, aPos, p1, p2);
    };

    computeVariant(k0al_alltrue_val, false, false, false, false, false, false);
    computeVariant(k0al_allreco_val, true, true, true, true, true, true);

    computeVariant(k0al_cvreco_val, true, false, false, false, false, false);
    computeVariant(k0al_avreco_val, false, true, false, false, false, false);
    computeVariant(k0al_vtxreco_val, true, true, false, false, false, false);

    computeVariant(k0al_d1magreco_val, false, false, true, false, false, false);
    computeVariant(k0al_d2magreco_val, false, false, false, false, true, false);
    computeVariant(k0al_d1dirreco_val, false, false, false, true, false, false);
    computeVariant(k0al_d2dirreco_val, false, false, false, false, false, true);
    computeVariant(k0al_d1preco_val, false, false, true, true, false, false);
    computeVariant(k0al_d2preco_val, false, false, false, false, true, true);

    computeVariant(k0al_d12magreco_val, false, false, true, false, true, false);
    computeVariant(k0al_d12dirreco_val, false, false, false, true, false, true);
    computeVariant(k0al_d12preco_val, false, false, true, true, true, true);

    computeVariant(k0al_allreco_cvtrue_val, false, true, true, true, true, true);
    computeVariant(k0al_allreco_avtrue_val, true, false, true, true, true, true);
    computeVariant(k0al_allreco_d1true_val, true, true, false, false, true, true);
    computeVariant(k0al_allreco_d2true_val, true, true, true, true, false, false);
  }

  output.FillVectorVar(k0al_alltrue, k0al_alltrue_val);
  output.FillVectorVar(k0al_allreco, k0al_allreco_val);
  output.FillVectorVar(k0al_cvreco, k0al_cvreco_val);
  output.FillVectorVar(k0al_avreco, k0al_avreco_val);
  output.FillVectorVar(k0al_vtxreco, k0al_vtxreco_val);
  output.FillVectorVar(k0al_d1magreco, k0al_d1magreco_val);
  output.FillVectorVar(k0al_d2magreco, k0al_d2magreco_val);
  output.FillVectorVar(k0al_d1dirreco, k0al_d1dirreco_val);
  output.FillVectorVar(k0al_d2dirreco, k0al_d2dirreco_val);
  output.FillVectorVar(k0al_d1preco, k0al_d1preco_val);
  output.FillVectorVar(k0al_d2preco, k0al_d2preco_val);
  output.FillVectorVar(k0al_d12magreco, k0al_d12magreco_val);
  output.FillVectorVar(k0al_d12dirreco, k0al_d12dirreco_val);
  output.FillVectorVar(k0al_d12preco, k0al_d12preco_val);
  output.FillVectorVar(k0al_allreco_cvtrue, k0al_allreco_cvtrue_val);
  output.FillVectorVar(k0al_allreco_avtrue, k0al_allreco_avtrue_val);
  output.FillVectorVar(k0al_allreco_d1true, k0al_allreco_d1true_val);
  output.FillVectorVar(k0al_allreco_d2true, k0al_allreco_d2true_val);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Parent(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
  Int_t k0partruepdg_val = -999;
  Float_t k0partrueendmom_val = -999.0f;
  Float_t k0partruelength_val = -999.0f;
  Float_t k0parrecolength_val = -999.0f;
  Int_t k0parndau_val = -1;
  Int_t k0parisbeam_val = -1;

  if (candidate) {
    AnaParticlePD* recoParent = nullptr;
    if (candidate->CreationVertex && candidate->CreationVertex->BeamParticle) {
      recoParent = static_cast<AnaParticlePD*>(candidate->CreationVertex->BeamParticle);
    } else if (candidate->Parent) {
      recoParent = static_cast<AnaParticlePD*>(candidate->Parent);
    }
    AnaTrueParticlePD* trueParent =
        (recoParent && recoParent->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoParent->TrueObject) : nullptr;

    if (trueParent) {
      k0partruepdg_val = trueParent->PDG;
      k0partrueendmom_val = trueParent->MomentumEnd;

      const TVector3 trueStart(trueParent->Position[0], trueParent->Position[1], trueParent->Position[2]);
      const TVector3 trueEnd(trueParent->PositionEnd[0], trueParent->PositionEnd[1], trueParent->PositionEnd[2]);
      if (IsValidPos(trueStart) && IsValidPos(trueEnd)) {
        k0partruelength_val = static_cast<Float_t>((trueEnd - trueStart).Mag());
      }
    }

    if (recoParent) {
      k0parrecolength_val = recoParent->Length;
      k0parndau_val = static_cast<Int_t>(recoParent->Daughters.size());
      k0parisbeam_val = recoParent->isPandora ? 1 : 0;
    }
  }

  output.FillVectorVar(k0partruepdg, k0partruepdg_val);
  output.FillVectorVar(k0partrueendmom, k0partrueendmom_val);
  output.FillVectorVar(k0partruelength, k0partruelength_val);
  output.FillVectorVar(k0parrecolength, k0parrecolength_val);
  output.FillVectorVar(k0parndau, k0parndau_val);
  output.FillVectorVar(k0parisbeam, k0parisbeam_val);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0CreationVtx(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
  Float_t k0cvtxpandoraresidual_val = -999.0f;
  Float_t k0cvtxfitresidual_val = -999.0f;
  Float_t k0cvtxpandoraresidualx_val = -999.0f;
  Float_t k0cvtxpandoraresidualy_val = -999.0f;
  Float_t k0cvtxpandoraresidualz_val = -999.0f;
  Float_t k0cvtxfitresidualx_val = -999.0f;
  Float_t k0cvtxfitresidualy_val = -999.0f;
  Float_t k0cvtxfitresidualz_val = -999.0f;
  Float_t k0protonmomentumreco_val = -999.0f;
  Float_t k0protonmomentumtrue_val = -999.0f;
  Int_t k0protontruepdg_val = -999;
  Float_t k0protonchi2ndfproton_val = -999.0f;
  Float_t k0protontruelength_val = -999.0f;
  Float_t k0protonrecolength_val = -999.0f;
  Int_t k0hasproton_val = 0;
  Int_t k0creationdegeneracy_val = -999;

  if (candidate) {
    AnaCreationVertexPD* creationVtx = candidate->CreationVertex;
    AnaParticlePD* recoParent = nullptr;
    if (creationVtx && creationVtx->BeamParticle) {
      recoParent = creationVtx->BeamParticle;
    } else if (candidate->Parent) {
      recoParent = static_cast<AnaParticlePD*>(candidate->Parent);
    }

    TVector3 trueCreation(-999.0, -999.0, -999.0);
    bool hasTrueCreation = false;
    AnaTrueParticlePD* trueK0 = candidate->TrueObject ? static_cast<AnaTrueParticlePD*>(candidate->TrueObject) : nullptr;
    if (trueK0) {
      trueCreation.SetXYZ(trueK0->Position[0], trueK0->Position[1], trueK0->Position[2]);
      hasTrueCreation = IsValidPos(trueCreation);
    }
    if (!hasTrueCreation) {
      AnaTrueParticlePD* trueParent =
          (recoParent && recoParent->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoParent->TrueObject) : nullptr;
      if (trueParent) {
        trueCreation.SetXYZ(trueParent->PositionEnd[0], trueParent->PositionEnd[1], trueParent->PositionEnd[2]);
        hasTrueCreation = IsValidPos(trueCreation);
      }
    }

    const bool hasRecoPandora = recoParent &&
                                recoParent->PositionEnd[0] > -900.f &&
                                recoParent->PositionEnd[1] > -900.f &&
                                recoParent->PositionEnd[2] > -900.f;
    const bool hasRecoFit = creationVtx &&
                            creationVtx->Position[0] > -900.f &&
                            creationVtx->Position[1] > -900.f &&
                            creationVtx->Position[2] > -900.f;

    if (hasTrueCreation && hasRecoPandora) {
      const TVector3 recoPandora(recoParent->PositionEnd[0], recoParent->PositionEnd[1], recoParent->PositionEnd[2]);
      k0cvtxpandoraresidual_val = static_cast<Float_t>((recoPandora - trueCreation).Mag());
      k0cvtxpandoraresidualx_val = recoParent->PositionEnd[0] - static_cast<Float_t>(trueCreation.X());
      k0cvtxpandoraresidualy_val = recoParent->PositionEnd[1] - static_cast<Float_t>(trueCreation.Y());
      k0cvtxpandoraresidualz_val = recoParent->PositionEnd[2] - static_cast<Float_t>(trueCreation.Z());
    }

    if (hasTrueCreation && hasRecoFit) {
      const TVector3 recoFit(creationVtx->Position[0], creationVtx->Position[1], creationVtx->Position[2]);
      k0cvtxfitresidual_val = static_cast<Float_t>((recoFit - trueCreation).Mag());
      k0cvtxfitresidualx_val = creationVtx->Position[0] - static_cast<Float_t>(trueCreation.X());
      k0cvtxfitresidualy_val = creationVtx->Position[1] - static_cast<Float_t>(trueCreation.Y());
      k0cvtxfitresidualz_val = creationVtx->Position[2] - static_cast<Float_t>(trueCreation.Z());
    }

    if (creationVtx && creationVtx->SecondParticle) {
      AnaParticlePD* secondParticle = creationVtx->SecondParticle;
      const bool secondIsValid =
          secondParticle->PositionStart[0] > -900.f &&
          secondParticle->PositionStart[1] > -900.f &&
          secondParticle->PositionStart[2] > -900.f;
      k0hasproton_val = secondIsValid ? 1 : 0;
      k0protonmomentumreco_val = secondParticle->Momentum;
      k0protonrecolength_val = secondParticle->Length;
      k0protonchi2ndfproton_val = pdAnaUtils::Chi2PIDChi2PerHit(secondParticle, 2212);

      AnaTrueParticlePD* secondTrue =
          secondParticle->TrueObject ? static_cast<AnaTrueParticlePD*>(secondParticle->TrueObject) : nullptr;
      if (secondTrue) {
        k0protonmomentumtrue_val = secondTrue->Momentum;
        k0protontruepdg_val = secondTrue->PDG;
        const TVector3 trueStart(secondTrue->Position[0], secondTrue->Position[1], secondTrue->Position[2]);
        const TVector3 trueEnd(secondTrue->PositionEnd[0], secondTrue->PositionEnd[1], secondTrue->PositionEnd[2]);
        if (IsValidPos(trueStart) && IsValidPos(trueEnd)) {
          k0protontruelength_val = static_cast<Float_t>((trueEnd - trueStart).Mag());
        }
      }
    }

    if (creationVtx) {
      k0creationdegeneracy_val = creationVtx->Degeneracy;
    }
  }

  output.FillVectorVar(k0cvtxpandoraresidual, k0cvtxpandoraresidual_val);
  output.FillVectorVar(k0cvtxfitresidual, k0cvtxfitresidual_val);
  output.FillVectorVar(k0cvtxpandoraresidualx, k0cvtxpandoraresidualx_val);
  output.FillVectorVar(k0cvtxpandoraresidualy, k0cvtxpandoraresidualy_val);
  output.FillVectorVar(k0cvtxpandoraresidualz, k0cvtxpandoraresidualz_val);
  output.FillVectorVar(k0cvtxfitresidualx, k0cvtxfitresidualx_val);
  output.FillVectorVar(k0cvtxfitresidualy, k0cvtxfitresidualy_val);
  output.FillVectorVar(k0cvtxfitresidualz, k0cvtxfitresidualz_val);
  output.FillVectorVar(k0protonmomentumreco, k0protonmomentumreco_val);
  output.FillVectorVar(k0protonmomentumtrue, k0protonmomentumtrue_val);
  output.FillVectorVar(k0protontruepdg, k0protontruepdg_val);
  output.FillVectorVar(k0protonchi2ndfproton, k0protonchi2ndfproton_val);
  output.FillVectorVar(k0protontruelength, k0protontruelength_val);
  output.FillVectorVar(k0protonrecolength, k0protonrecolength_val);
  output.FillVectorVar(k0hasproton, k0hasproton_val);
  output.FillVectorVar(k0creationdegeneracy, k0creationdegeneracy_val);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Kinematics(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
  // Initialize all variables to sentinel
  Float_t k0truecreationx_val = -999.0f;
  Float_t k0truecreationy_val = -999.0f;
  Float_t k0truecreationz_val = -999.0f;
  Float_t k0trueannihilationx_val = -999.0f;
  Float_t k0trueannihilationy_val = -999.0f;
  Float_t k0trueannihilationz_val = -999.0f;
  Float_t k0truemomentum_val = -999.0f;
  Float_t k0truedirectionx_val = -999.0f;
  Float_t k0truedirectiony_val = -999.0f;
  Float_t k0truedirectionz_val = -999.0f;

  Float_t k0creationpandorax_val = -999.0f;
  Float_t k0creationpandoray_val = -999.0f;
  Float_t k0creationpandoraz_val = -999.0f;

  Float_t k0creationfitx_val = -999.0f;
  Float_t k0creationfity_val = -999.0f;
  Float_t k0creationfitz_val = -999.0f;

  Float_t k0directionpandorax_val = -999.0f;
  Float_t k0directionpandoray_val = -999.0f;
  Float_t k0directionpandoraz_val = -999.0f;

  Float_t k0directionfitx_val = -999.0f;
  Float_t k0directionfity_val = -999.0f;
  Float_t k0directionfitz_val = -999.0f;

  Float_t k0vtxdau1truedirectionx_val = -999.0f;
  Float_t k0vtxdau1truedirectiony_val = -999.0f;
  Float_t k0vtxdau1truedirectionz_val = -999.0f;
  Float_t k0vtxdau2truedirectionx_val = -999.0f;
  Float_t k0vtxdau2truedirectiony_val = -999.0f;
  Float_t k0vtxdau2truedirectionz_val = -999.0f;

  Float_t k0vtxdau1directionpandorax_val = -999.0f;
  Float_t k0vtxdau1directionpandoray_val = -999.0f;
  Float_t k0vtxdau1directionpandoraz_val = -999.0f;
  Float_t k0vtxdau2directionpandorax_val = -999.0f;
  Float_t k0vtxdau2directionpandoray_val = -999.0f;
  Float_t k0vtxdau2directionpandoraz_val = -999.0f;

  Float_t k0vtxdau1directionfitx_val = -999.0f;
  Float_t k0vtxdau1directionfity_val = -999.0f;
  Float_t k0vtxdau1directionfitz_val = -999.0f;
  Float_t k0vtxdau2directionfitx_val = -999.0f;
  Float_t k0vtxdau2directionfity_val = -999.0f;
  Float_t k0vtxdau2directionfitz_val = -999.0f;

  if (candidate) {
    // Get reco parent
    AnaCreationVertexPD* creationVtx = candidate->CreationVertex;
    AnaParticlePD* recoParent = nullptr;
    if (creationVtx && creationVtx->BeamParticle) {
      recoParent = creationVtx->BeamParticle;
    } else if (candidate->Parent) {
      recoParent = static_cast<AnaParticlePD*>(candidate->Parent);
    }

    AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;

    // Get true K0
    AnaTrueParticlePD* trueK0 = candidate->TrueObject ? static_cast<AnaTrueParticlePD*>(candidate->TrueObject) : nullptr;
    AnaTrueParticlePD* trueParent = (recoParent && recoParent->TrueObject) ? static_cast<AnaTrueParticlePD*>(recoParent->TrueObject) : nullptr;

    // ===== TRUE K0 KINEMATICS =====
    TVector3 trueCreation(-999.0, -999.0, -999.0);
    bool hasTrueCreation = false;
    if (trueK0) {
      trueCreation.SetXYZ(trueK0->Position[0], trueK0->Position[1], trueK0->Position[2]);
      hasTrueCreation = IsValidPos(trueCreation);
    }
    if (!hasTrueCreation && trueParent) {
      trueCreation.SetXYZ(trueParent->PositionEnd[0], trueParent->PositionEnd[1], trueParent->PositionEnd[2]);
      hasTrueCreation = IsValidPos(trueCreation);
    }

    TVector3 trueAnnihilation(-999.0, -999.0, -999.0);
    bool hasTrueAnnihilation = false;
    if (trueK0) {
      trueAnnihilation.SetXYZ(trueK0->PositionEnd[0], trueK0->PositionEnd[1], trueK0->PositionEnd[2]);
      hasTrueAnnihilation = IsValidPos(trueAnnihilation);
    }

    if (hasTrueCreation) {
      k0truecreationx_val = static_cast<Float_t>(trueCreation.X());
      k0truecreationy_val = static_cast<Float_t>(trueCreation.Y());
      k0truecreationz_val = static_cast<Float_t>(trueCreation.Z());
    }

    if (hasTrueAnnihilation) {
      k0trueannihilationx_val = static_cast<Float_t>(trueAnnihilation.X());
      k0trueannihilationy_val = static_cast<Float_t>(trueAnnihilation.Y());
      k0trueannihilationz_val = static_cast<Float_t>(trueAnnihilation.Z());
    }

    if (trueK0 && trueK0->Momentum > 0.f) {
      k0truemomentum_val = trueK0->Momentum;
    }

    if (hasTrueCreation && hasTrueAnnihilation) {
      TVector3 trueDir = trueAnnihilation - trueCreation;
      if (trueDir.Mag() > 1e-6) {
        trueDir = trueDir.Unit();
        k0truedirectionx_val = static_cast<Float_t>(trueDir.X());
        k0truedirectiony_val = static_cast<Float_t>(trueDir.Y());
        k0truedirectionz_val = static_cast<Float_t>(trueDir.Z());
      }
    }

    // ===== RECO K0 CREATION POSITIONS =====
    if (recoParent) {
      if (recoParent->PositionEnd[0] > -900.f && recoParent->PositionEnd[1] > -900.f && recoParent->PositionEnd[2] > -900.f) {
        k0creationpandorax_val = recoParent->PositionEnd[0];
        k0creationpandoray_val = recoParent->PositionEnd[1];
        k0creationpandoraz_val = recoParent->PositionEnd[2];
      }
    }

    if (creationVtx) {
      if (creationVtx->Position[0] > -900.f && creationVtx->Position[1] > -900.f && creationVtx->Position[2] > -900.f) {
        k0creationfitx_val = creationVtx->Position[0];
        k0creationfity_val = creationVtx->Position[1];
        k0creationfitz_val = creationVtx->Position[2];
      }
    }

    // ===== RECO K0 DIRECTIONS =====
    if (vertex && recoParent) {
      TVector3 recoPandoraCreation(recoParent->PositionEnd[0], recoParent->PositionEnd[1], recoParent->PositionEnd[2]);
      TVector3 recoPandoraAnnihilation(vertex->PositionPandora[0], vertex->PositionPandora[1], vertex->PositionPandora[2]);
      if (IsValidPos(recoPandoraCreation) && IsValidPos(recoPandoraAnnihilation)) {
        TVector3 pandoraDir = recoPandoraAnnihilation - recoPandoraCreation;
        if (pandoraDir.Mag() > 1e-6) {
          pandoraDir = pandoraDir.Unit();
          k0directionpandorax_val = static_cast<Float_t>(pandoraDir.X());
          k0directionpandoray_val = static_cast<Float_t>(pandoraDir.Y());
          k0directionpandoraz_val = static_cast<Float_t>(pandoraDir.Z());
        }
      }
    }

    if (vertex && creationVtx) {
      TVector3 recoFitCreation(creationVtx->Position[0], creationVtx->Position[1], creationVtx->Position[2]);
      TVector3 recoFitAnnihilation(vertex->PositionFit[0], vertex->PositionFit[1], vertex->PositionFit[2]);
      if (IsValidPos(recoFitCreation) && IsValidPos(recoFitAnnihilation)) {
        TVector3 fitDir = recoFitAnnihilation - recoFitCreation;
        if (fitDir.Mag() > 1e-6) {
          fitDir = fitDir.Unit();
          k0directionfitx_val = static_cast<Float_t>(fitDir.X());
          k0directionfity_val = static_cast<Float_t>(fitDir.Y());
          k0directionfitz_val = static_cast<Float_t>(fitDir.Z());
        }
      }
    }

    // ===== TRUE DAUGHTER DIRECTIONS =====
    if (vertex && vertex->Particles.size() > 0) {
      AnaParticlePD* recoDau1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      AnaTrueParticlePD* trueDau1 = recoDau1 ? static_cast<AnaTrueParticlePD*>(recoDau1->TrueObject) : nullptr;
      if (trueDau1 && trueDau1->Direction[0] > -900.f && trueDau1->Direction[1] > -900.f && trueDau1->Direction[2] > -900.f) {
        k0vtxdau1truedirectionx_val = trueDau1->Direction[0];
        k0vtxdau1truedirectiony_val = trueDau1->Direction[1];
        k0vtxdau1truedirectionz_val = trueDau1->Direction[2];
      }
    }

    if (vertex && vertex->Particles.size() > 1) {
      AnaParticlePD* recoDau2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      AnaTrueParticlePD* trueDau2 = recoDau2 ? static_cast<AnaTrueParticlePD*>(recoDau2->TrueObject) : nullptr;
      if (trueDau2 && trueDau2->Direction[0] > -900.f && trueDau2->Direction[1] > -900.f && trueDau2->Direction[2] > -900.f) {
        k0vtxdau2truedirectionx_val = trueDau2->Direction[0];
        k0vtxdau2truedirectiony_val = trueDau2->Direction[1];
        k0vtxdau2truedirectionz_val = trueDau2->Direction[2];
      }
    }

    // ===== RECO DAUGHTER DIRECTIONS: PANDORA =====
    if (vertex && vertex->Particles.size() > 0) {
      AnaParticlePD* recoDau1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      if (recoDau1 && recoDau1->DirectionStart[0] > -900.f && recoDau1->DirectionStart[1] > -900.f && recoDau1->DirectionStart[2] > -900.f) {
        k0vtxdau1directionpandorax_val = recoDau1->DirectionStart[0];
        k0vtxdau1directionpandoray_val = recoDau1->DirectionStart[1];
        k0vtxdau1directionpandoraz_val = recoDau1->DirectionStart[2];
      }
    }

    if (vertex && vertex->Particles.size() > 1) {
      AnaParticlePD* recoDau2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      if (recoDau2 && recoDau2->DirectionStart[0] > -900.f && recoDau2->DirectionStart[1] > -900.f && recoDau2->DirectionStart[2] > -900.f) {
        k0vtxdau2directionpandorax_val = recoDau2->DirectionStart[0];
        k0vtxdau2directionpandoray_val = recoDau2->DirectionStart[1];
        k0vtxdau2directionpandoraz_val = recoDau2->DirectionStart[2];
      }
    }

    // ===== RECO DAUGHTER DIRECTIONS: FIT-EXTRAPOLATED =====
    const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
    const double trackFitDistanceFromStart = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
    if (vertex && vertex->Particles.size() > 0) {
      AnaParticlePD* recoDau1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      if (recoDau1) {
        std::vector<double> fit1;
        pdAnaUtils::ExtrapolateTrack(recoDau1, fit1, trackFitLength, true, trackFitDistanceFromStart);
        const bool fit1Valid = (fit1.size() >= 6 && fit1[3] > -900.0 && fit1[4] > -900.0 && fit1[5] > -900.0 &&
                                std::isfinite(fit1[3]) && std::isfinite(fit1[4]) && std::isfinite(fit1[5]));
        if (fit1Valid) {
          k0vtxdau1directionfitx_val = static_cast<Float_t>(fit1[3]);
          k0vtxdau1directionfity_val = static_cast<Float_t>(fit1[4]);
          k0vtxdau1directionfitz_val = static_cast<Float_t>(fit1[5]);
        }
      }
    }

    if (vertex && vertex->Particles.size() > 1) {
      AnaParticlePD* recoDau2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      if (recoDau2) {
        std::vector<double> fit2;
        pdAnaUtils::ExtrapolateTrack(recoDau2, fit2, trackFitLength, true, trackFitDistanceFromStart);
        const bool fit2Valid = (fit2.size() >= 6 && fit2[3] > -900.0 && fit2[4] > -900.0 && fit2[5] > -900.0 &&
                                std::isfinite(fit2[3]) && std::isfinite(fit2[4]) && std::isfinite(fit2[5]));
        if (fit2Valid) {
          k0vtxdau2directionfitx_val = static_cast<Float_t>(fit2[3]);
          k0vtxdau2directionfity_val = static_cast<Float_t>(fit2[4]);
          k0vtxdau2directionfitz_val = static_cast<Float_t>(fit2[5]);
        }
      }
    }
  }

  // Fill all variables
  output.FillVectorVar(k0truecreationx, k0truecreationx_val);
  output.FillVectorVar(k0truecreationy, k0truecreationy_val);
  output.FillVectorVar(k0truecreationz, k0truecreationz_val);
  output.FillVectorVar(k0trueannihilationx, k0trueannihilationx_val);
  output.FillVectorVar(k0trueannihilationy, k0trueannihilationy_val);
  output.FillVectorVar(k0trueannihilationz, k0trueannihilationz_val);
  output.FillVectorVar(k0truemomentum, k0truemomentum_val);
  output.FillVectorVar(k0truedirectionx, k0truedirectionx_val);
  output.FillVectorVar(k0truedirectiony, k0truedirectiony_val);
  output.FillVectorVar(k0truedirectionz, k0truedirectionz_val);

  output.FillVectorVar(k0creationpandorax, k0creationpandorax_val);
  output.FillVectorVar(k0creationpandoray, k0creationpandoray_val);
  output.FillVectorVar(k0creationpandoraz, k0creationpandoraz_val);

  output.FillVectorVar(k0creationfitx, k0creationfitx_val);
  output.FillVectorVar(k0creationfity, k0creationfity_val);
  output.FillVectorVar(k0creationfitz, k0creationfitz_val);

  output.FillVectorVar(k0directionpandorax, k0directionpandorax_val);
  output.FillVectorVar(k0directionpandoray, k0directionpandoray_val);
  output.FillVectorVar(k0directionpandoraz, k0directionpandoraz_val);

  output.FillVectorVar(k0directionfitx, k0directionfitx_val);
  output.FillVectorVar(k0directionfity, k0directionfity_val);
  output.FillVectorVar(k0directionfitz, k0directionfitz_val);

  output.FillVectorVar(k0vtxdau1truedirectionx, k0vtxdau1truedirectionx_val);
  output.FillVectorVar(k0vtxdau1truedirectiony, k0vtxdau1truedirectiony_val);
  output.FillVectorVar(k0vtxdau1truedirectionz, k0vtxdau1truedirectionz_val);
  output.FillVectorVar(k0vtxdau2truedirectionx, k0vtxdau2truedirectionx_val);
  output.FillVectorVar(k0vtxdau2truedirectiony, k0vtxdau2truedirectiony_val);
  output.FillVectorVar(k0vtxdau2truedirectionz, k0vtxdau2truedirectionz_val);

  output.FillVectorVar(k0vtxdau1directionpandorax, k0vtxdau1directionpandorax_val);
  output.FillVectorVar(k0vtxdau1directionpandoray, k0vtxdau1directionpandoray_val);
  output.FillVectorVar(k0vtxdau1directionpandoraz, k0vtxdau1directionpandoraz_val);
  output.FillVectorVar(k0vtxdau2directionpandorax, k0vtxdau2directionpandorax_val);
  output.FillVectorVar(k0vtxdau2directionpandoray, k0vtxdau2directionpandoray_val);
  output.FillVectorVar(k0vtxdau2directionpandoraz, k0vtxdau2directionpandoraz_val);

  output.FillVectorVar(k0vtxdau1directionfitx, k0vtxdau1directionfitx_val);
  output.FillVectorVar(k0vtxdau1directionfity, k0vtxdau1directionfity_val);
  output.FillVectorVar(k0vtxdau1directionfitz, k0vtxdau1directionfitz_val);
  output.FillVectorVar(k0vtxdau2directionfitx, k0vtxdau2directionfitx_val);
  output.FillVectorVar(k0vtxdau2directionfity, k0vtxdau2directionfity_val);
  output.FillVectorVar(k0vtxdau2directionfitz, k0vtxdau2directionfitz_val);
}

//********************************************************************
//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaAnnihilationVertexPD* vertex,
                                                     const AnaEventB& event, Int_t excludedParentUniqueID,
                                                     AnaCreationVertexPD* creationVertex){
    // Fill all variables for a single K0 vertex
  Float_t invalidPos[3] = {-999.0f, -999.0f, -999.0f};
  Float_t k0vtxtruepos_val[3] = {-999.0f, -999.0f, -999.0f};
  Float_t k0vtxtrueoriginaldistance_val = -999.0f;
  Float_t k0vtxdegeneracy_val = -999.0f;
  Float_t k0vtxoriginaldistance_val = -999.0f;
  Float_t k0vtxpandoraresidual_val = -999.0f;
  Float_t k0vtxfitresidual_val = -999.0f;
  Float_t k0vtxpandorax_val = -999.0f;
  Float_t k0vtxpandoray_val = -999.0f;
  Float_t k0vtxpandoraz_val = -999.0f;
  Float_t k0vtxfitx_val = -999.0f;
  Float_t k0vtxfity_val = -999.0f;
  Float_t k0vtxfitz_val = -999.0f;
  Float_t k0vtxpandoraresidualx_val = -999.0f;
  Float_t k0vtxpandoraresidualy_val = -999.0f;
  Float_t k0vtxpandoraresidualz_val = -999.0f;
  Float_t k0vtxfitresidualx_val = -999.0f;
  Float_t k0vtxfitresidualy_val = -999.0f;
  Float_t k0vtxfitresidualz_val = -999.0f;
  Float_t k0vtxmomentum_val = -999.0f;
  Float_t k0vtxinvariantmass_val = -999.0f;
  Float_t k0vtxmomentumpandora_val = -999.0f;
  Float_t k0vtxinvariantmasspandora_val = -999.0f;
  Float_t k0vtxopeninganglepandora_val = -999.0f;
  Float_t k0vtxmomentumfit_val = -999.0f;
  Float_t k0vtxinvariantmassfit_val = -999.0f;
  Float_t k0vtxopeninganglefit_val = -999.0f;
  Float_t k0vtxresultantmomentumreco_val = -999.0f;
  Float_t k0vtxresultantmomentumtrue_val = -999.0f;
  Int_t k0vtxjointk0smomused_val = -1;
  Float_t k0vtxjointsigmap1_val = -999.0f;
  Float_t k0vtxjointsigmap2_val = -999.0f;
  Float_t k0vtxjointsigmam_val = -999.0f;
  Float_t k0vtxjointdmdp1_val = -999.0f;
  Float_t k0vtxjointdmdp2_val = -999.0f;
  Float_t k0vtxjointrconstraint_val = -999.0f;
  Float_t k0vtxjointdeltachidedx_val = -999.0f;
  Int_t k0vtxjointdebugclass_val = 0;

  if(vertex){
    k0vtxoriginaldistance_val = vertex->OriginalDistance;
    k0vtxdegeneracy_val = static_cast<Float_t>(
        pdAnnihilationUtils::ComputeAnnihilationVertexDegeneracyWithExclusion(
            event, vertex, excludedParentUniqueID, creationVertex));
    if (vertex->PositionPandora[0] > -900.f && vertex->PositionPandora[1] > -900.f && vertex->PositionPandora[2] > -900.f) {
      k0vtxpandorax_val = vertex->PositionPandora[0];
      k0vtxpandoray_val = vertex->PositionPandora[1];
      k0vtxpandoraz_val = vertex->PositionPandora[2];
    }
    if (vertex->PositionFit[0] > -900.f && vertex->PositionFit[1] > -900.f && vertex->PositionFit[2] > -900.f) {
      k0vtxfitx_val = vertex->PositionFit[0];
      k0vtxfity_val = vertex->PositionFit[1];
      k0vtxfitz_val = vertex->PositionFit[2];
    }
    k0vtxmomentum_val = vertex->Momentum;
    k0vtxinvariantmass_val = vertex->InvariantMass;
    k0vtxmomentumpandora_val = vertex->MomentumPandora;
    k0vtxinvariantmasspandora_val = vertex->InvariantMassPandora;
    k0vtxmomentumfit_val = vertex->MomentumFit;
    k0vtxinvariantmassfit_val = vertex->InvariantMassFit;
    k0vtxresultantmomentumreco_val = vertex->Momentum;
    k0vtxjointk0smomused_val = vertex->JointK0sMomentumUsed;
    k0vtxjointsigmap1_val = vertex->JointK0sSigmaP1GeV;
    k0vtxjointsigmap2_val = vertex->JointK0sSigmaP2GeV;
    k0vtxjointsigmam_val = vertex->JointK0sSigmaMEventGeV;
    k0vtxjointdmdp1_val = vertex->JointK0sDmDp1;
    k0vtxjointdmdp2_val = vertex->JointK0sDmDp2;
    k0vtxjointrconstraint_val = vertex->JointK0sMomentumConstraintRatioR;
    k0vtxjointdeltachidedx_val = vertex->JointK0sMomentumDedxChi2Degradation;
    k0vtxjointdebugclass_val = vertex->JointK0sDebugClass;

    if(vertex->Particles.size() >= 2) {
      AnaParticlePD* recoParticle1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      AnaParticlePD* recoParticle2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      AnaTrueParticlePD* trueParticle1 = recoParticle1 ? static_cast<AnaTrueParticlePD*>(recoParticle1->TrueObject) : nullptr;
      AnaTrueParticlePD* trueParticle2 = recoParticle2 ? static_cast<AnaTrueParticlePD*>(recoParticle2->TrueObject) : nullptr;

      if (recoParticle1 && recoParticle2) {
        const TVector3 dirPandora1(recoParticle1->DirectionStart[0], recoParticle1->DirectionStart[1], recoParticle1->DirectionStart[2]);
        const TVector3 dirPandora2(recoParticle2->DirectionStart[0], recoParticle2->DirectionStart[1], recoParticle2->DirectionStart[2]);
        k0vtxopeninganglepandora_val = AngleBetweenDirections(dirPandora1, dirPandora2);

        TVector3 dirFit1 = dirPandora1;
        TVector3 dirFit2 = dirPandora2;
        const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
        const double trackFitDistanceFromStart =
            ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
        std::vector<double> fit1;
        std::vector<double> fit2;
        pdAnaUtils::ExtrapolateTrack(recoParticle1, fit1, trackFitLength, true, trackFitDistanceFromStart);
        pdAnaUtils::ExtrapolateTrack(recoParticle2, fit2, trackFitLength, true, trackFitDistanceFromStart);
        const bool fit1Valid = (fit1.size() >= 6 && fit1[3] > -900.0 && fit1[4] > -900.0 && fit1[5] > -900.0 &&
                                std::isfinite(fit1[3]) && std::isfinite(fit1[4]) && std::isfinite(fit1[5]));
        const bool fit2Valid = (fit2.size() >= 6 && fit2[3] > -900.0 && fit2[4] > -900.0 && fit2[5] > -900.0 &&
                                std::isfinite(fit2[3]) && std::isfinite(fit2[4]) && std::isfinite(fit2[5]));
        if (fit1Valid) dirFit1.SetXYZ(fit1[3], fit1[4], fit1[5]);
        if (fit2Valid) dirFit2.SetXYZ(fit2[3], fit2[4], fit2[5]);
        k0vtxopeninganglefit_val = AngleBetweenDirections(dirFit1, dirFit2);
      }

      if(trueParticle1 && trueParticle2) {
        const Float_t dx = trueParticle1->Position[0] - trueParticle2->Position[0];
        const Float_t dy = trueParticle1->Position[1] - trueParticle2->Position[1];
        const Float_t dz = trueParticle1->Position[2] - trueParticle2->Position[2];
        k0vtxtrueoriginaldistance_val = sqrt(dx*dx + dy*dy + dz*dz);
        TVector3 dirTrue1(trueParticle1->Direction[0], trueParticle1->Direction[1], trueParticle1->Direction[2]);
        TVector3 dirTrue2(trueParticle2->Direction[0], trueParticle2->Direction[1], trueParticle2->Direction[2]);
        if (dirTrue1.Mag2() > 0.0 && dirTrue2.Mag2() > 0.0 && trueParticle1->Momentum > 0.0f && trueParticle2->Momentum > 0.0f) {
          const TVector3 pTrueTot = trueParticle1->Momentum * dirTrue1.Unit() + trueParticle2->Momentum * dirTrue2.Unit();
          k0vtxresultantmomentumtrue_val = static_cast<Float_t>(pTrueTot.Mag());
        }

        // True vertex is valid only when both true daughters are consistent with a common start point.
        constexpr Float_t kSameTrueVertexTolerance = 1e-2f;
        if (k0vtxtrueoriginaldistance_val <= kSameTrueVertexTolerance) {
          k0vtxtruepos_val[0] = 0.5f * static_cast<Float_t>(trueParticle1->Position[0] + trueParticle2->Position[0]);
          k0vtxtruepos_val[1] = 0.5f * static_cast<Float_t>(trueParticle1->Position[1] + trueParticle2->Position[1]);
          k0vtxtruepos_val[2] = 0.5f * static_cast<Float_t>(trueParticle1->Position[2] + trueParticle2->Position[2]);
        }
      }
    }

    if (k0vtxtruepos_val[0] > -900.f && k0vtxtruepos_val[1] > -900.f && k0vtxtruepos_val[2] > -900.f) {
      const TVector3 truePos(k0vtxtruepos_val[0], k0vtxtruepos_val[1], k0vtxtruepos_val[2]);

      if (vertex->PositionPandora[0] > -900.f && vertex->PositionPandora[1] > -900.f && vertex->PositionPandora[2] > -900.f) {
        const TVector3 pandoraPos(vertex->PositionPandora[0], vertex->PositionPandora[1], vertex->PositionPandora[2]);
        k0vtxpandoraresidual_val = static_cast<Float_t>((pandoraPos - truePos).Mag());
        k0vtxpandoraresidualx_val = vertex->PositionPandora[0] - k0vtxtruepos_val[0];
        k0vtxpandoraresidualy_val = vertex->PositionPandora[1] - k0vtxtruepos_val[1];
        k0vtxpandoraresidualz_val = vertex->PositionPandora[2] - k0vtxtruepos_val[2];
      }

      if (vertex->PositionFit[0] > -900.f && vertex->PositionFit[1] > -900.f && vertex->PositionFit[2] > -900.f) {
        const TVector3 fitPos(vertex->PositionFit[0], vertex->PositionFit[1], vertex->PositionFit[2]);
        k0vtxfitresidual_val = static_cast<Float_t>((fitPos - truePos).Mag());
        k0vtxfitresidualx_val = vertex->PositionFit[0] - k0vtxtruepos_val[0];
        k0vtxfitresidualy_val = vertex->PositionFit[1] - k0vtxtruepos_val[1];
        k0vtxfitresidualz_val = vertex->PositionFit[2] - k0vtxtruepos_val[2];
      }
    }

    output.FillMatrixVarFromArray(k0vtxpandorapos, vertex->PositionPandora, 3);
    output.FillMatrixVarFromArray(k0vtxfitpos, vertex->PositionFit, 3);

  } else {
    output.FillMatrixVarFromArray(k0vtxpandorapos, invalidPos, 3);
    output.FillMatrixVarFromArray(k0vtxfitpos, invalidPos, 3);
  }

  output.FillMatrixVarFromArray(k0vtxtruepos, k0vtxtruepos_val, 3);
  output.FillVectorVar(k0vtxoriginaldistance, k0vtxoriginaldistance_val);
  output.FillVectorVar(k0vtxtrueoriginaldistance, k0vtxtrueoriginaldistance_val);
  output.FillVectorVar(k0vtxdegeneracy, k0vtxdegeneracy_val);
  output.FillVectorVar(k0vtxpandorax, k0vtxpandorax_val);
  output.FillVectorVar(k0vtxpandoray, k0vtxpandoray_val);
  output.FillVectorVar(k0vtxpandoraz, k0vtxpandoraz_val);
  output.FillVectorVar(k0vtxfitx, k0vtxfitx_val);
  output.FillVectorVar(k0vtxfity, k0vtxfity_val);
  output.FillVectorVar(k0vtxfitz, k0vtxfitz_val);
  output.FillVectorVar(k0vtxpandoraresidual, k0vtxpandoraresidual_val);
  output.FillVectorVar(k0vtxfitresidual, k0vtxfitresidual_val);
  output.FillVectorVar(k0vtxpandoraresidualx, k0vtxpandoraresidualx_val);
  output.FillVectorVar(k0vtxpandoraresidualy, k0vtxpandoraresidualy_val);
  output.FillVectorVar(k0vtxpandoraresidualz, k0vtxpandoraresidualz_val);
  output.FillVectorVar(k0vtxfitresidualx, k0vtxfitresidualx_val);
  output.FillVectorVar(k0vtxfitresidualy, k0vtxfitresidualy_val);
  output.FillVectorVar(k0vtxfitresidualz, k0vtxfitresidualz_val);
  output.FillVectorVar(k0vtxmomentum, k0vtxmomentum_val);
  output.FillVectorVar(k0vtxinvariantmass, k0vtxinvariantmass_val);
  output.FillVectorVar(k0vtxmomentumpandora, k0vtxmomentumpandora_val);
  output.FillVectorVar(k0vtxinvariantmasspandora, k0vtxinvariantmasspandora_val);
  output.FillVectorVar(k0vtxopeninganglepandora, k0vtxopeninganglepandora_val);
  output.FillVectorVar(k0vtxmomentumfit, k0vtxmomentumfit_val);
  output.FillVectorVar(k0vtxinvariantmassfit, k0vtxinvariantmassfit_val);
  output.FillVectorVar(k0vtxopeninganglefit, k0vtxopeninganglefit_val);
  output.FillVectorVar(k0vtxresultantmomentumreco, k0vtxresultantmomentumreco_val);
  output.FillVectorVar(k0vtxresultantmomentumtrue, k0vtxresultantmomentumtrue_val);
  output.FillVectorVar(k0vtxjointk0smomused, k0vtxjointk0smomused_val);
  output.FillVectorVar(k0vtxjointsigmap1, k0vtxjointsigmap1_val);
  output.FillVectorVar(k0vtxjointsigmap2, k0vtxjointsigmap2_val);
  output.FillVectorVar(k0vtxjointsigmam, k0vtxjointsigmam_val);
  output.FillVectorVar(k0vtxjointdmdp1, k0vtxjointdmdp1_val);
  output.FillVectorVar(k0vtxjointdmdp2, k0vtxjointdmdp2_val);
  output.FillVectorVar(k0vtxjointrconstraint, k0vtxjointrconstraint_val);
  output.FillVectorVar(k0vtxjointdeltachidedx, k0vtxjointdeltachidedx_val);
  output.FillVectorVar(k0vtxjointdebugclass, k0vtxjointdebugclass_val);
  FillNeutralKaonVariables_K0vtxDaughters(output, vertex, event);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtxDaughters(OutputManager& output, AnaAnnihilationVertexPD* vertex, const AnaEventB& event){
//********************************************************************
  Float_t k0vtxdau1momentumreco_val = -999.0f;
  Float_t k0vtxdau2momentumreco_val = -999.0f;
  Int_t k0vtxdau1mommethod_val = -1;
  Int_t k0vtxdau2mommethod_val = -1;
  Float_t k0vtxdau1extchi2ndf_val = -999.0f;
  Float_t k0vtxdau2extchi2ndf_val = -999.0f;
  Float_t k0vtxdau1dedxdrift_val = -999.0f;
  Float_t k0vtxdau2dedxdrift_val = -999.0f;
  Float_t k0vtxdau1dedxsigma_val = -999.0f;
  Float_t k0vtxdau2dedxsigma_val = -999.0f;
  Int_t k0vtxdau1dedxfitok_val = -1;
  Int_t k0vtxdau2dedxfitok_val = -1;
  Float_t k0vtxdau1momentumtrue_val = -999.0f;
  Float_t k0vtxdau2momentumtrue_val = -999.0f;
  Int_t k0vtxdau1trueendproc_val = -1;
  Int_t k0vtxdau2trueendproc_val = -1;
  Float_t k0vtxdau1trueendmom_val = -999.0f;
  Float_t k0vtxdau2trueendmom_val = -999.0f;
  Float_t k0vtxdau1truestartmom_val = -999.0f;
  Float_t k0vtxdau2truestartmom_val = -999.0f;
  Float_t k0vtxdau1truestartx_val = -999.0f;
  Float_t k0vtxdau1truestarty_val = -999.0f;
  Float_t k0vtxdau1truestartz_val = -999.0f;
  Float_t k0vtxdau2truestartx_val = -999.0f;
  Float_t k0vtxdau2truestarty_val = -999.0f;
  Float_t k0vtxdau2truestartz_val = -999.0f;
  Int_t k0vtxdau1truendau_val = -1;
  Int_t k0vtxdau2truendau_val = -1;
  Float_t k0vtxdau1truelength_val = -999.0f;
  Float_t k0vtxdau2truelength_val = -999.0f;
  Float_t k0vtxdau1recolength_val = -999.0f;
  Float_t k0vtxdau2recolength_val = -999.0f;
  Float_t k0vtxdau1recostartx_val = -999.0f;
  Float_t k0vtxdau1recostarty_val = -999.0f;
  Float_t k0vtxdau1recostartz_val = -999.0f;
  Float_t k0vtxdau2recostartx_val = -999.0f;
  Float_t k0vtxdau2recostarty_val = -999.0f;
  Float_t k0vtxdau2recostartz_val = -999.0f;
  Int_t k0vtxdau1nhitsreco_val = -1;
  Int_t k0vtxdau2nhitsreco_val = -1;
  Float_t k0vtxdau1protonchi2ndf_val = -999.0f;
  Float_t k0vtxdau2protonchi2ndf_val = -999.0f;
  Float_t k0vtxdau1pionchi2ndf_val = -999.0f;
  Float_t k0vtxdau2pionchi2ndf_val = -999.0f;
  Float_t k0vtxdau1kaonchi2ndf_val = -999.0f;
  Float_t k0vtxdau2kaonchi2ndf_val = -999.0f;
  Int_t k0vtxdau1truepdg_val = -999;
  Int_t k0vtxdau2truepdg_val = -999;
  Int_t k0vtxdau1ndaughtersreco_val = -1;
  Int_t k0vtxdau2ndaughtersreco_val = -1;
  Int_t k0vtxdau1nrecodau_val = -1;
  Int_t k0vtxdau2nrecodau_val = -1;
  Float_t k0vtxdau1recovisiblee_val = -999.0f;
  Float_t k0vtxdau2recovisiblee_val = -999.0f;
  Float_t k0vtxdaughtersrecovisiblee_val = -999.0f;

  if (vertex) {
    k0vtxdau1mommethod_val = vertex->Daughter1MomentumMethod;
    k0vtxdau2mommethod_val = vertex->Daughter2MomentumMethod;
    k0vtxdau1extchi2ndf_val = vertex->Daughter1ExtensionChi2Ndf;
    k0vtxdau2extchi2ndf_val = vertex->Daughter2ExtensionChi2Ndf;
    k0vtxdau1dedxdrift_val = vertex->Daughter1ExtensionDedxBias;
    k0vtxdau2dedxdrift_val = vertex->Daughter2ExtensionDedxBias;
    k0vtxdau1dedxsigma_val = vertex->Daughter1ExtensionDedxSigma;
    k0vtxdau2dedxsigma_val = vertex->Daughter2ExtensionDedxSigma;
    k0vtxdau1dedxfitok_val = vertex->Daughter1ExtensionDedxFitOk;
    k0vtxdau2dedxfitok_val = vertex->Daughter2ExtensionDedxFitOk;

    if (vertex->Particles.size() >= 2) {
      AnaParticlePD* recoParticle1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      AnaParticlePD* recoParticle2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      AnaTrueParticlePD* trueParticle1 = recoParticle1 ? static_cast<AnaTrueParticlePD*>(recoParticle1->TrueObject) : nullptr;
      AnaTrueParticlePD* trueParticle2 = recoParticle2 ? static_cast<AnaTrueParticlePD*>(recoParticle2->TrueObject) : nullptr;
      if (recoParticle1) {
        k0vtxdau1momentumreco_val = recoParticle1->Momentum;
        k0vtxdau1recolength_val = recoParticle1->Length;
        k0vtxdau1nhitsreco_val = recoParticle1->NHitsPerPlane[2];
        const TVector3 recoStart1(recoParticle1->PositionStart[0], recoParticle1->PositionStart[1], recoParticle1->PositionStart[2]);
        if (IsValidPos(recoStart1)) {
          k0vtxdau1recostartx_val = static_cast<Float_t>(recoStart1.X());
          k0vtxdau1recostarty_val = static_cast<Float_t>(recoStart1.Y());
          k0vtxdau1recostartz_val = static_cast<Float_t>(recoStart1.Z());
        }
        k0vtxdau1protonchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle1, 2212);
        k0vtxdau1pionchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle1, 211);
        k0vtxdau1kaonchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle1, 321);
        k0vtxdau1ndaughtersreco_val = static_cast<Int_t>(recoParticle1->Daughters.size());
      }
      if (recoParticle2) {
        k0vtxdau2momentumreco_val = recoParticle2->Momentum;
        k0vtxdau2recolength_val = recoParticle2->Length;
        k0vtxdau2nhitsreco_val = recoParticle2->NHitsPerPlane[2];
        const TVector3 recoStart2(recoParticle2->PositionStart[0], recoParticle2->PositionStart[1], recoParticle2->PositionStart[2]);
        if (IsValidPos(recoStart2)) {
          k0vtxdau2recostartx_val = static_cast<Float_t>(recoStart2.X());
          k0vtxdau2recostarty_val = static_cast<Float_t>(recoStart2.Y());
          k0vtxdau2recostartz_val = static_cast<Float_t>(recoStart2.Z());
        }
        k0vtxdau2protonchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle2, 2212);
        k0vtxdau2pionchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle2, 211);
        k0vtxdau2kaonchi2ndf_val = pdAnaUtils::Chi2PIDChi2PerHit(recoParticle2, 321);
        k0vtxdau2ndaughtersreco_val = static_cast<Int_t>(recoParticle2->Daughters.size());
      }
      if (trueParticle1) {
        k0vtxdau1truepdg_val = trueParticle1->PDG;
        k0vtxdau1momentumtrue_val = trueParticle1->Momentum;
        k0vtxdau1truestartmom_val = trueParticle1->Momentum;
        k0vtxdau1trueendmom_val = trueParticle1->MomentumEnd;
        k0vtxdau1trueendproc_val = static_cast<Int_t>(trueParticle1->ProcessEnd);
        k0vtxdau1truendau_val = static_cast<Int_t>(trueParticle1->Daughters.size());
        k0vtxdau1nrecodau_val = CountRecoObjectsForTrueDaughters(trueParticle1, event);
        const TVector3 trueStart1(trueParticle1->Position[0], trueParticle1->Position[1], trueParticle1->Position[2]);
        const TVector3 trueEnd1(trueParticle1->PositionEnd[0], trueParticle1->PositionEnd[1], trueParticle1->PositionEnd[2]);
        if (IsValidPos(trueStart1)) {
          k0vtxdau1truestartx_val = static_cast<Float_t>(trueStart1.X());
          k0vtxdau1truestarty_val = static_cast<Float_t>(trueStart1.Y());
          k0vtxdau1truestartz_val = static_cast<Float_t>(trueStart1.Z());
        }
        if (IsValidPos(trueStart1) && IsValidPos(trueEnd1)) {
          k0vtxdau1truelength_val = static_cast<Float_t>((trueEnd1 - trueStart1).Mag());
        }
      }
      if (trueParticle2) {
        k0vtxdau2truepdg_val = trueParticle2->PDG;
        k0vtxdau2momentumtrue_val = trueParticle2->Momentum;
        k0vtxdau2truestartmom_val = trueParticle2->Momentum;
        k0vtxdau2trueendmom_val = trueParticle2->MomentumEnd;
        k0vtxdau2trueendproc_val = static_cast<Int_t>(trueParticle2->ProcessEnd);
        k0vtxdau2truendau_val = static_cast<Int_t>(trueParticle2->Daughters.size());
        k0vtxdau2nrecodau_val = CountRecoObjectsForTrueDaughters(trueParticle2, event);
        const TVector3 trueStart2(trueParticle2->Position[0], trueParticle2->Position[1], trueParticle2->Position[2]);
        const TVector3 trueEnd2(trueParticle2->PositionEnd[0], trueParticle2->PositionEnd[1], trueParticle2->PositionEnd[2]);
        if (IsValidPos(trueStart2)) {
          k0vtxdau2truestartx_val = static_cast<Float_t>(trueStart2.X());
          k0vtxdau2truestarty_val = static_cast<Float_t>(trueStart2.Y());
          k0vtxdau2truestartz_val = static_cast<Float_t>(trueStart2.Z());
        }
        if (IsValidPos(trueStart2) && IsValidPos(trueEnd2)) {
          k0vtxdau2truelength_val = static_cast<Float_t>((trueEnd2 - trueStart2).Mag());
        }
      }

      double vis1MeV = 0.;
      double vis2MeV = 0.;
      if (recoParticle1) {
        vis1MeV = SumRecoDaughterSubtreesVisibleEnergyMeV(recoParticle1);
        k0vtxdau1recovisiblee_val = static_cast<Float_t>(vis1MeV / 1000.0);
      }
      if (recoParticle2) {
        vis2MeV = SumRecoDaughterSubtreesVisibleEnergyMeV(recoParticle2);
        k0vtxdau2recovisiblee_val = static_cast<Float_t>(vis2MeV / 1000.0);
      }
      k0vtxdaughtersrecovisiblee_val = static_cast<Float_t>((vis1MeV + vis2MeV) / 1000.0);
    }
  }

  output.FillVectorVar(k0vtxdau1momentumreco, k0vtxdau1momentumreco_val);
  output.FillVectorVar(k0vtxdau2momentumreco, k0vtxdau2momentumreco_val);
  output.FillVectorVar(k0vtxdau1mommethod, k0vtxdau1mommethod_val);
  output.FillVectorVar(k0vtxdau2mommethod, k0vtxdau2mommethod_val);
  output.FillVectorVar(k0vtxdau1extchi2ndf, k0vtxdau1extchi2ndf_val);
  output.FillVectorVar(k0vtxdau2extchi2ndf, k0vtxdau2extchi2ndf_val);
  output.FillVectorVar(k0vtxdau1dedxdrift, k0vtxdau1dedxdrift_val);
  output.FillVectorVar(k0vtxdau2dedxdrift, k0vtxdau2dedxdrift_val);
  output.FillVectorVar(k0vtxdau1dedxsigma, k0vtxdau1dedxsigma_val);
  output.FillVectorVar(k0vtxdau2dedxsigma, k0vtxdau2dedxsigma_val);
  output.FillVectorVar(k0vtxdau1dedxfitok, k0vtxdau1dedxfitok_val);
  output.FillVectorVar(k0vtxdau2dedxfitok, k0vtxdau2dedxfitok_val);
  output.FillVectorVar(k0vtxdau1momentumtrue, k0vtxdau1momentumtrue_val);
  output.FillVectorVar(k0vtxdau2momentumtrue, k0vtxdau2momentumtrue_val);
  output.FillVectorVar(k0vtxdau1trueendproc, k0vtxdau1trueendproc_val);
  output.FillVectorVar(k0vtxdau2trueendproc, k0vtxdau2trueendproc_val);
  output.FillVectorVar(k0vtxdau1trueendmom, k0vtxdau1trueendmom_val);
  output.FillVectorVar(k0vtxdau2trueendmom, k0vtxdau2trueendmom_val);
  output.FillVectorVar(k0vtxdau1truestartmom, k0vtxdau1truestartmom_val);
  output.FillVectorVar(k0vtxdau2truestartmom, k0vtxdau2truestartmom_val);
  output.FillVectorVar(k0vtxdau1truestartx, k0vtxdau1truestartx_val);
  output.FillVectorVar(k0vtxdau1truestarty, k0vtxdau1truestarty_val);
  output.FillVectorVar(k0vtxdau1truestartz, k0vtxdau1truestartz_val);
  output.FillVectorVar(k0vtxdau2truestartx, k0vtxdau2truestartx_val);
  output.FillVectorVar(k0vtxdau2truestarty, k0vtxdau2truestarty_val);
  output.FillVectorVar(k0vtxdau2truestartz, k0vtxdau2truestartz_val);
  output.FillVectorVar(k0vtxdau1truendau, k0vtxdau1truendau_val);
  output.FillVectorVar(k0vtxdau2truendau, k0vtxdau2truendau_val);
  output.FillVectorVar(k0vtxdau1truelength, k0vtxdau1truelength_val);
  output.FillVectorVar(k0vtxdau2truelength, k0vtxdau2truelength_val);
  output.FillVectorVar(k0vtxdau1recolength, k0vtxdau1recolength_val);
  output.FillVectorVar(k0vtxdau2recolength, k0vtxdau2recolength_val);
  output.FillVectorVar(k0vtxdau1recostartx, k0vtxdau1recostartx_val);
  output.FillVectorVar(k0vtxdau1recostarty, k0vtxdau1recostarty_val);
  output.FillVectorVar(k0vtxdau1recostartz, k0vtxdau1recostartz_val);
  output.FillVectorVar(k0vtxdau2recostartx, k0vtxdau2recostartx_val);
  output.FillVectorVar(k0vtxdau2recostarty, k0vtxdau2recostarty_val);
  output.FillVectorVar(k0vtxdau2recostartz, k0vtxdau2recostartz_val);
  output.FillVectorVar(k0vtxdau1nhitsreco, k0vtxdau1nhitsreco_val);
  output.FillVectorVar(k0vtxdau2nhitsreco, k0vtxdau2nhitsreco_val);
  output.FillVectorVar(k0vtxdau1protonchi2ndf, k0vtxdau1protonchi2ndf_val);
  output.FillVectorVar(k0vtxdau2protonchi2ndf, k0vtxdau2protonchi2ndf_val);
  output.FillVectorVar(k0vtxdau1pionchi2ndf, k0vtxdau1pionchi2ndf_val);
  output.FillVectorVar(k0vtxdau2pionchi2ndf, k0vtxdau2pionchi2ndf_val);
  output.FillVectorVar(k0vtxdau1kaonchi2ndf, k0vtxdau1kaonchi2ndf_val);
  output.FillVectorVar(k0vtxdau2kaonchi2ndf, k0vtxdau2kaonchi2ndf_val);
  output.FillVectorVar(k0vtxdau1truepdg, k0vtxdau1truepdg_val);
  output.FillVectorVar(k0vtxdau2truepdg, k0vtxdau2truepdg_val);
  output.FillVectorVar(k0vtxdau1ndaughtersreco, k0vtxdau1ndaughtersreco_val);
  output.FillVectorVar(k0vtxdau2ndaughtersreco, k0vtxdau2ndaughtersreco_val);
  output.FillVectorVar(k0vtxdau1nrecodau, k0vtxdau1nrecodau_val);
  output.FillVectorVar(k0vtxdau2nrecodau, k0vtxdau2nrecodau_val);
  output.FillVectorVar(k0vtxdau1recovisiblee, k0vtxdau1recovisiblee_val);
  output.FillVectorVar(k0vtxdau2recovisiblee, k0vtxdau2recovisiblee_val);
  output.FillVectorVar(k0vtxdaughtersrecovisiblee, k0vtxdaughtersrecovisiblee_val);
}

//********************************************************************
