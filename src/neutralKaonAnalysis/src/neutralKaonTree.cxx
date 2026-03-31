#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdDataClasses.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include "AnalysisUtils.hxx"
#include "ParticleId.hxx"
#include "HEPConstants.hxx"
#include "TH2F.h"
#include "TProfile.h"
#include <algorithm>
#include <cmath>

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

  bool IsValidPos(const TVector3& p) {
    return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) &&
           p.X() > -900.f && p.Y() > -900.f && p.Z() > -900.f;
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
}

//********************************************************************
void neutralKaonTree::WriteHitDistanceProfiles(OutputManager& output){
//********************************************************************
  TTree* tree = output.GetTree();
  if (!tree) {
    return;
  }
  TFile* file = tree->GetCurrentFile();
  if (!file) {
    return;
  }

  file->cd();
  for (Int_t i = 0; i < 2; ++i) {
    if (gK0DauHitDistVsTravel2D[i]) {
      gK0DauHitDistVsTravel2D[i]->Write("", TObject::kOverwrite);
    }
    if (gK0DauHitDistVsTravelProf[i]) {
      gK0DauHitDistVsTravelProf[i]->Write("", TObject::kOverwrite);
    }
  }
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Particle(OutputManager& output, UInt_t nmax){
//********************************************************************
  AddVarMaxSizeVF(output, k0lengthpandora, "Neutral length using annihilation Pandora position", nk0, nmax);
  AddVarMaxSizeVF(output, k0lengthfit, "Neutral length using annihilation Fit position", nk0, nmax);
  AddVarMaxSizeVF(output, k0alignmentpandora,
                   "K0 alignment (Pandora): angle (rad) between creation→annihilation(Pandora) and vertex Σp (Pandora dirs)", nk0,
                   nmax);
  AddVarMaxSizeVF(output, k0alignmentfit,
                   "K0 alignment (Fit): angle (rad) between creation→annihilation(Fit) and vertex Σp (fit dirs)", nk0,
                   nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax){

  AddVarMaxSizeVI(output, k0nvtxbeforefiltering, "Number of annihilation vertices before overlap filtering", nk0, nmax);
  AddVarMaxSizeVI(output, k0nvtxafterfiltering, "Number of annihilation vertices after overlap filtering", nk0, nmax);

  //Vertex system variables
  AddVarMaxSize3MF(output, k0vtxtruepos, "K0 vertex true position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxoriginaldistance, "K0 vertex original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueoriginaldistance, "K0 vertex true original distance", nk0, nmax);
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
  AddVarMaxSizeVF(output, k0vtxmomentumfit, "K0 vertex momentum using fit directions", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxinvariantmassfit, "K0 vertex invariant mass using fit directions", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                               Int_t nVerticesBeforeFiltering, Int_t nVerticesAfterFiltering,
                                               AnaBeamB* beam){
  (void)beam;
  output.FillVectorVar(k0nvtxbeforefiltering, nVerticesBeforeFiltering);
  output.FillVectorVar(k0nvtxafterfiltering, nVerticesAfterFiltering);
  neutralKaonTree::FillNeutralKaonVariables_K0Particle(output, candidate);

  if(candidate){
    AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
    neutralKaonTree::FillNeutralKaonVariables_K0vtx(output, vertex);

    if (vertex && neutralKaonAnaUtils::IsSignalCandidate(candidate, event)) {
      if (vertex->Particles.size() > 0 && vertex->Particles[0]) {
        AnaParticlePD* recoDau1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
        AnaTrueParticlePD* trueDau1 = recoDau1 ? static_cast<AnaTrueParticlePD*>(recoDau1->TrueObject) : nullptr;
        FillHitDistanceToTrueLineProfiles(output, recoDau1, trueDau1, 0);
      }
      if (vertex->Particles.size() > 1 && vertex->Particles[1]) {
        AnaParticlePD* recoDau2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
        AnaTrueParticlePD* trueDau2 = recoDau2 ? static_cast<AnaTrueParticlePD*>(recoDau2->TrueObject) : nullptr;
        FillHitDistanceToTrueLineProfiles(output, recoDau2, trueDau2, 1);
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
  Float_t alignmentPandora = -999.0f;
  Float_t alignmentFit = -999.0f;

  if (candidate) {
    lengthPandora = candidate->LengthPandora;
    lengthFit = candidate->LengthFit;
    alignmentPandora = candidate->AlignmentPandora;
    alignmentFit = candidate->AlignmentFit;
  }

  output.FillVectorVar(k0lengthpandora, lengthPandora);
  output.FillVectorVar(k0lengthfit, lengthFit);
  output.FillVectorVar(k0alignmentpandora, alignmentPandora);
  output.FillVectorVar(k0alignmentfit, alignmentFit);
}

//********************************************************************
//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaAnnihilationVertexPD* vertex){
    // Fill all variables for a single K0 vertex
  Float_t invalidPos[3] = {-999.0f, -999.0f, -999.0f};
  Float_t k0vtxtruepos_val[3] = {-999.0f, -999.0f, -999.0f};
  Float_t k0vtxtrueoriginaldistance_val = -999.0f;
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
  Float_t k0vtxmomentumfit_val = -999.0f;
  Float_t k0vtxinvariantmassfit_val = -999.0f;

  if(vertex){
    k0vtxoriginaldistance_val = vertex->OriginalDistance;
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

    if(vertex->Particles.size() >= 2) {
      AnaParticlePD* recoParticle1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
      AnaParticlePD* recoParticle2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
      AnaTrueParticlePD* trueParticle1 = recoParticle1 ? static_cast<AnaTrueParticlePD*>(recoParticle1->TrueObject) : nullptr;
      AnaTrueParticlePD* trueParticle2 = recoParticle2 ? static_cast<AnaTrueParticlePD*>(recoParticle2->TrueObject) : nullptr;

      if(trueParticle1 && trueParticle2) {
        const Float_t dx = trueParticle1->Position[0] - trueParticle2->Position[0];
        const Float_t dy = trueParticle1->Position[1] - trueParticle2->Position[1];
        const Float_t dz = trueParticle1->Position[2] - trueParticle2->Position[2];
        k0vtxtrueoriginaldistance_val = sqrt(dx*dx + dy*dy + dz*dz);

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
  output.FillVectorVar(k0vtxmomentumfit, k0vtxmomentumfit_val);
  output.FillVectorVar(k0vtxinvariantmassfit, k0vtxinvariantmassfit_val);
}

//********************************************************************
