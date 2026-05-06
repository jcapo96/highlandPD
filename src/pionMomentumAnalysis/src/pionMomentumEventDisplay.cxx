#include "pionMomentumEventDisplay.hxx"

#include "OutputManager.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "pdDataClasses.hxx"
#include "pdMomReconstructionMCS.hxx"

#include <TEveElement.h>
#include <TEveLine.h>
#include <TEvePointSet.h>
#include <TEveScene.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

/// Same `trajTruthForTrueMcs` resolution as `pionMomentumAnalysis.cxx` (true MCS uses this object's TrjPoints).
const AnaTrueParticlePD* TrajTruthForTrueMcs(const AnaParticlePD* mainTrack, const AnaParticlePD* beamParticle) {
  if (!mainTrack) return nullptr;
  const AnaTrueParticlePD* mainTrueParticle =
      mainTrack->TrueObject ? static_cast<const AnaTrueParticlePD*>(mainTrack->TrueObject) : nullptr;
  const AnaTrueParticlePD* beamTrueParticle =
      (beamParticle && beamParticle->TrueObject)
          ? static_cast<const AnaTrueParticlePD*>(beamParticle->TrueObject)
          : nullptr;
  const AnaTrueParticlePD* trajTruthForTrueMcs = nullptr;
  if (mainTrueParticle && beamTrueParticle && mainTrueParticle->ID == beamTrueParticle->ID &&
      mainTrueParticle->TrjPoints.empty() && !beamTrueParticle->TrjPoints.empty()) {
    trajTruthForTrueMcs = beamTrueParticle;
  } else if (mainTrueParticle) {
    trajTruthForTrueMcs = mainTrueParticle;
  } else {
    trajTruthForTrueMcs = beamTrueParticle;
  }
  if (trajTruthForTrueMcs && trajTruthForTrueMcs->TrjPoints.empty() && beamTrueParticle &&
      !beamTrueParticle->TrjPoints.empty()) {
    trajTruthForTrueMcs = beamTrueParticle;
  }
  return trajTruthForTrueMcs;
}

} // namespace

ClassImp(pionMomentumEventDisplay);

const Int_t pionMomentumEventDisplay::kMaxTrjPoints;
const Int_t pionMomentumEventDisplay::kMaxMcsSegments;

//********************************************************************
pionMomentumEventDisplay::pionMomentumEventDisplay() : pdEventDisplay() {
//********************************************************************
  _eventDisplayClassName = "pionMomentumEventDisplay";
  _nTrjPoints = 0;
  _trj_X = new Float_t[kMaxTrjPoints];
  _trj_Y = new Float_t[kMaxTrjPoints];
  _trj_Z = new Float_t[kMaxTrjPoints];
  _trj_dEdx = new Float_t[kMaxTrjPoints];

  _nMcsSegments = 0;
  _mcs_mainTrack_uniqueID = -1;
  _edHasMcsBranches = false;
  _nMcsTrueSegments = 0;
  _mcs_true_traj_trueID = -999;
  _edHasMcsTrueBranches = false;
  for (Int_t s = 0; s < kMaxMcsSegments; ++s) {
    _mcs_seg_chord_X1[s] = -999.f; _mcs_seg_chord_Y1[s] = -999.f; _mcs_seg_chord_Z1[s] = -999.f;
    _mcs_seg_chord_X2[s] = -999.f; _mcs_seg_chord_Y2[s] = -999.f; _mcs_seg_chord_Z2[s] = -999.f;
    _mcs_seg_fit_X1[s]   = -999.f; _mcs_seg_fit_Y1[s]   = -999.f; _mcs_seg_fit_Z1[s]   = -999.f;
    _mcs_seg_fit_X2[s]   = -999.f; _mcs_seg_fit_Y2[s]   = -999.f; _mcs_seg_fit_Z2[s]   = -999.f;
    _mcs_true_seg_chord_X1[s] = -999.f; _mcs_true_seg_chord_Y1[s] = -999.f; _mcs_true_seg_chord_Z1[s] = -999.f;
    _mcs_true_seg_chord_X2[s] = -999.f; _mcs_true_seg_chord_Y2[s] = -999.f; _mcs_true_seg_chord_Z2[s] = -999.f;
    _mcs_true_seg_fit_X1[s]   = -999.f; _mcs_true_seg_fit_Y1[s]   = -999.f; _mcs_true_seg_fit_Z1[s]   = -999.f;
    _mcs_true_seg_fit_X2[s]   = -999.f; _mcs_true_seg_fit_Y2[s]   = -999.f; _mcs_true_seg_fit_Z2[s]   = -999.f;
  }
}

//********************************************************************
pionMomentumEventDisplay::~pionMomentumEventDisplay() {
//********************************************************************
  delete[] _trj_X;
  delete[] _trj_Y;
  delete[] _trj_Z;
  delete[] _trj_dEdx;
}

//********************************************************************
void pionMomentumEventDisplay::AddAnalysisVariables(OutputManager& output, Int_t tree_index) {
//********************************************************************
  output.AddVectorVar(tree_index, edparticle_nTrjPoints, "ED_particle_nTrjPoints", "I",
                      "Number of trajectory points per particle", ednParticles, "ED_nParticles", -kMaxParticles);
  output.AddVectorVar(tree_index, edparticle_trjPointIndex, "ED_particle_trjPointIndex", "I",
                      "Starting trajectory-point index for each particle", ednParticles, "ED_nParticles", -kMaxParticles);

  output.AddVectorVar(tree_index, edtrj_X, "ED_trj_X", "F", "Trajectory point X positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_Y, "ED_trj_Y", "F", "Trajectory point Y positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_Z, "ED_trj_Z", "F", "Trajectory point Z positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_dEdx, "ED_trj_dEdx", "F", "Trajectory point dEdx (placeholder)",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);

  output.AddVar(tree_index, edmcs_mainTrack_uniqueID, "ED_mcs_mainTrack_uniqueID", "I",
                "UniqueID of the reco particle whose MCS segments are stored");
  output.AddVectorVar(tree_index, edmcs_seg_chord_X1, "ED_mcs_seg_chord_X1", "F",
                      "MCS segment chord start X (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_chord_Y1, "ED_mcs_seg_chord_Y1", "F",
                      "MCS segment chord start Y (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_chord_Z1, "ED_mcs_seg_chord_Z1", "F",
                      "MCS segment chord start Z (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_chord_X2, "ED_mcs_seg_chord_X2", "F",
                      "MCS segment chord end X (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_chord_Y2, "ED_mcs_seg_chord_Y2", "F",
                      "MCS segment chord end Y (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_chord_Z2, "ED_mcs_seg_chord_Z2", "F",
                      "MCS segment chord end Z (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_X1, "ED_mcs_seg_fit_X1", "F",
                      "MCS segment PCA-fit line start X (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_Y1, "ED_mcs_seg_fit_Y1", "F",
                      "MCS segment PCA-fit line start Y (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_Z1, "ED_mcs_seg_fit_Z1", "F",
                      "MCS segment PCA-fit line start Z (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_X2, "ED_mcs_seg_fit_X2", "F",
                      "MCS segment PCA-fit line end X (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_Y2, "ED_mcs_seg_fit_Y2", "F",
                      "MCS segment PCA-fit line end Y (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_seg_fit_Z2, "ED_mcs_seg_fit_Z2", "F",
                      "MCS segment PCA-fit line end Z (reco-SCE)", ednMcsSegments, "ED_nMcsSegments", -kMaxMcsSegments);

  output.AddVar(tree_index, edmcs_true_traj_trueID, "ED_mcs_true_traj_trueID", "I",
                "True particle ID whose trajectory was used for true MCS (SCE-distorted G4 points)");
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_X1, "ED_mcs_true_seg_chord_X1", "F",
                      "True MCS segment chord start X (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_Y1, "ED_mcs_true_seg_chord_Y1", "F",
                      "True MCS segment chord start Y (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_Z1, "ED_mcs_true_seg_chord_Z1", "F",
                      "True MCS segment chord start Z (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_X2, "ED_mcs_true_seg_chord_X2", "F",
                      "True MCS segment chord end X (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_Y2, "ED_mcs_true_seg_chord_Y2", "F",
                      "True MCS segment chord end Y (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_chord_Z2, "ED_mcs_true_seg_chord_Z2", "F",
                      "True MCS segment chord end Z (SCE-distorted G4)", ednMcsTrueSegments, "ED_nMcsTrueSegments",
                      -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_X1, "ED_mcs_true_seg_fit_X1", "F",
                      "True MCS segment PCA-fit line start X (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_Y1, "ED_mcs_true_seg_fit_Y1", "F",
                      "True MCS segment PCA-fit line start Y (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_Z1, "ED_mcs_true_seg_fit_Z1", "F",
                      "True MCS segment PCA-fit line start Z (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_X2, "ED_mcs_true_seg_fit_X2", "F",
                      "True MCS segment PCA-fit line end X (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_Y2, "ED_mcs_true_seg_fit_Y2", "F",
                      "True MCS segment PCA-fit line end Y (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
  output.AddVectorVar(tree_index, edmcs_true_seg_fit_Z2, "ED_mcs_true_seg_fit_Z2", "F",
                      "True MCS segment PCA-fit line end Z (SCE-distorted G4)", ednMcsTrueSegments,
                      "ED_nMcsTrueSegments", -kMaxMcsSegments);
}

//********************************************************************
void pionMomentumEventDisplay::FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
  Int_t trjIndex = 0;
  Int_t particleIndex = 0;
  for (Int_t i = 0; i < event.nParticles && i < kMaxParticles; ++i) {
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[i]);
    if (!particle) continue;

    const Int_t nTrj = static_cast<Int_t>(particle->TrjPoints.size());
    output.FillVectorVarForceIndex(edparticle_nTrjPoints, nTrj, particleIndex);
    output.FillVectorVarForceIndex(edparticle_trjPointIndex, trjIndex, particleIndex);

    for (Int_t t = 0; t < nTrj && trjIndex < kMaxTrjPoints; ++t) {
      const AnaTrajectoryPointPD& trj = particle->TrjPoints[nTrj - 1 - t];
      output.FillVectorVar(edtrj_X, static_cast<Float_t>(trj.Position.X()));
      output.FillVectorVar(edtrj_Y, static_cast<Float_t>(trj.Position.Y()));
      output.FillVectorVar(edtrj_Z, static_cast<Float_t>(trj.Position.Z()));
      output.FillVectorVar(edtrj_dEdx, static_cast<Float_t>(-999.0));
      output.IncrementCounter(ednTrjPoints);
      ++trjIndex;
    }
    ++particleIndex;
  }

  Int_t mcsMainUID = -1;
  Int_t mcsTrueTrajId = -999;
  const ToyBoxPD* tbox = static_cast<const ToyBoxPD*>(box);
  if (tbox && tbox->MainTrack) {
    mcsMainUID = tbox->MainTrack->UniqueID;

    pdMomReconstruction::MCSLikelihoodConfig cfg;
    if (ND::params().HasParameter("pionMomentumAnalysis.MCSRadiationLengthCm"))
      cfg.radiationLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSRadiationLengthCm");
    if (ND::params().HasParameter("pionMomentumAnalysis.MCSTargetSegmentLengthCm"))
      cfg.targetSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSTargetSegmentLengthCm");
    if (ND::params().HasParameter("pionMomentumAnalysis.MCSMinSegmentLengthCm"))
      cfg.minSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSMinSegmentLengthCm");
    if (ND::params().HasParameter("pionMomentumAnalysis.MCStheta0FloorRad"))
      cfg.theta0FloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCStheta0FloorRad");
    if (ND::params().HasParameter("pionMomentumAnalysis.MCSMaxAbsDeltaThetaRad"))
      cfg.maxAbsDeltaThetaRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSMaxAbsDeltaThetaRad");

    std::vector<pdMomReconstruction::MCSSegmentGeometry> segs;
    pdMomReconstruction::BuildPionMcsSegmentGeometry(*tbox->MainTrack, cfg, segs);

    Int_t nSeg = 0;
    for (size_t s = 0; s < segs.size() && nSeg < kMaxMcsSegments; ++s) {
      const auto& g = segs[s];
      const Float_t x1 = static_cast<Float_t>(g.startPoint.X());
      const Float_t y1 = static_cast<Float_t>(g.startPoint.Y());
      const Float_t z1 = static_cast<Float_t>(g.startPoint.Z());
      const Float_t x2 = static_cast<Float_t>(g.endPoint.X());
      const Float_t y2 = static_cast<Float_t>(g.endPoint.Y());
      const Float_t z2 = static_cast<Float_t>(g.endPoint.Z());
      const TVector3 dir = (g.fittedDirection.Mag2() > 1e-20) ? g.fittedDirection.Unit()
                                                              : (g.endPoint - g.startPoint).Unit();
      const double half = 0.5 * std::max(0.0, g.arcLengthCm);
      const TVector3 p1 = g.centroid - half * dir;
      const TVector3 p2 = g.centroid + half * dir;

      output.FillVectorVar(edmcs_seg_chord_X1, x1);
      output.FillVectorVar(edmcs_seg_chord_Y1, y1);
      output.FillVectorVar(edmcs_seg_chord_Z1, z1);
      output.FillVectorVar(edmcs_seg_chord_X2, x2);
      output.FillVectorVar(edmcs_seg_chord_Y2, y2);
      output.FillVectorVar(edmcs_seg_chord_Z2, z2);
      output.FillVectorVar(edmcs_seg_fit_X1, static_cast<Float_t>(p1.X()));
      output.FillVectorVar(edmcs_seg_fit_Y1, static_cast<Float_t>(p1.Y()));
      output.FillVectorVar(edmcs_seg_fit_Z1, static_cast<Float_t>(p1.Z()));
      output.FillVectorVar(edmcs_seg_fit_X2, static_cast<Float_t>(p2.X()));
      output.FillVectorVar(edmcs_seg_fit_Y2, static_cast<Float_t>(p2.Y()));
      output.FillVectorVar(edmcs_seg_fit_Z2, static_cast<Float_t>(p2.Z()));
      output.IncrementCounter(ednMcsSegments);
      ++nSeg;
    }

    AnaBeamPD* beam = event.Beam ? static_cast<AnaBeamPD*>(event.Beam) : nullptr;
    AnaParticlePD* beampart =
        (beam && beam->BeamParticle) ? static_cast<AnaParticlePD*>(beam->BeamParticle) : nullptr;
    const AnaTrueParticlePD* trajTruth = TrajTruthForTrueMcs(tbox->MainTrack, beampart);
    if (trajTruth) mcsTrueTrajId = trajTruth->ID;

    std::vector<TVector3> trueOrdered;
    if (trajTruth) {
      trueOrdered.reserve(trajTruth->TrjPoints.size());
      for (const auto& tp : trajTruth->TrjPoints) {
        if (!tp.IsInTPC) continue;
        const TVector3& p = tp.Position;
        if (std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900.0 &&
            p.Y() > -900.0 && p.Z() > -900.0)
          trueOrdered.push_back(p);
      }
    }
    std::vector<pdMomReconstruction::MCSSegmentGeometry> trueSegs;
    if (!trueOrdered.empty())
      pdMomReconstruction::BuildPionMcsSegmentGeometryFromOrderedPositions(trueOrdered, cfg, trueSegs);

    Int_t nTrueSeg = 0;
    for (size_t s = 0; s < trueSegs.size() && nTrueSeg < kMaxMcsSegments; ++s) {
      const auto& g = trueSegs[s];
      const Float_t x1 = static_cast<Float_t>(g.startPoint.X());
      const Float_t y1 = static_cast<Float_t>(g.startPoint.Y());
      const Float_t z1 = static_cast<Float_t>(g.startPoint.Z());
      const Float_t x2 = static_cast<Float_t>(g.endPoint.X());
      const Float_t y2 = static_cast<Float_t>(g.endPoint.Y());
      const Float_t z2 = static_cast<Float_t>(g.endPoint.Z());
      const TVector3 dir = (g.fittedDirection.Mag2() > 1e-20) ? g.fittedDirection.Unit()
                                                              : (g.endPoint - g.startPoint).Unit();
      const double half = 0.5 * std::max(0.0, g.arcLengthCm);
      const TVector3 p1 = g.centroid - half * dir;
      const TVector3 p2 = g.centroid + half * dir;

      output.FillVectorVar(edmcs_true_seg_chord_X1, x1);
      output.FillVectorVar(edmcs_true_seg_chord_Y1, y1);
      output.FillVectorVar(edmcs_true_seg_chord_Z1, z1);
      output.FillVectorVar(edmcs_true_seg_chord_X2, x2);
      output.FillVectorVar(edmcs_true_seg_chord_Y2, y2);
      output.FillVectorVar(edmcs_true_seg_chord_Z2, z2);
      output.FillVectorVar(edmcs_true_seg_fit_X1, static_cast<Float_t>(p1.X()));
      output.FillVectorVar(edmcs_true_seg_fit_Y1, static_cast<Float_t>(p1.Y()));
      output.FillVectorVar(edmcs_true_seg_fit_Z1, static_cast<Float_t>(p1.Z()));
      output.FillVectorVar(edmcs_true_seg_fit_X2, static_cast<Float_t>(p2.X()));
      output.FillVectorVar(edmcs_true_seg_fit_Y2, static_cast<Float_t>(p2.Y()));
      output.FillVectorVar(edmcs_true_seg_fit_Z2, static_cast<Float_t>(p2.Z()));
      output.IncrementCounter(ednMcsTrueSegments);
      ++nTrueSeg;
    }
  }
  output.FillVar(edmcs_mainTrack_uniqueID, mcsMainUID);
  output.FillVar(edmcs_true_traj_trueID, mcsTrueTrajId);
}

//********************************************************************
bool pionMomentumEventDisplay::ReadAnalysisData(TTree* tree) {
//********************************************************************
  tree->SetBranchAddress("ED_nTrjPoints", &_nTrjPoints);
  tree->SetBranchAddress("ED_particle_nTrjPoints", _particle_nTrjPoints);
  tree->SetBranchAddress("ED_particle_trjPointIndex", _particle_trjPointIndex);
  tree->SetBranchAddress("ED_trj_X", _trj_X);
  tree->SetBranchAddress("ED_trj_Y", _trj_Y);
  tree->SetBranchAddress("ED_trj_Z", _trj_Z);
  tree->SetBranchAddress("ED_trj_dEdx", _trj_dEdx);

  _edHasMcsBranches = tree->GetBranch("ED_nMcsSegments") != nullptr;
  if (_edHasMcsBranches) {
    tree->SetBranchAddress("ED_nMcsSegments", &_nMcsSegments);
    tree->SetBranchAddress("ED_mcs_mainTrack_uniqueID", &_mcs_mainTrack_uniqueID);
    tree->SetBranchAddress("ED_mcs_seg_chord_X1", _mcs_seg_chord_X1);
    tree->SetBranchAddress("ED_mcs_seg_chord_Y1", _mcs_seg_chord_Y1);
    tree->SetBranchAddress("ED_mcs_seg_chord_Z1", _mcs_seg_chord_Z1);
    tree->SetBranchAddress("ED_mcs_seg_chord_X2", _mcs_seg_chord_X2);
    tree->SetBranchAddress("ED_mcs_seg_chord_Y2", _mcs_seg_chord_Y2);
    tree->SetBranchAddress("ED_mcs_seg_chord_Z2", _mcs_seg_chord_Z2);
    tree->SetBranchAddress("ED_mcs_seg_fit_X1", _mcs_seg_fit_X1);
    tree->SetBranchAddress("ED_mcs_seg_fit_Y1", _mcs_seg_fit_Y1);
    tree->SetBranchAddress("ED_mcs_seg_fit_Z1", _mcs_seg_fit_Z1);
    tree->SetBranchAddress("ED_mcs_seg_fit_X2", _mcs_seg_fit_X2);
    tree->SetBranchAddress("ED_mcs_seg_fit_Y2", _mcs_seg_fit_Y2);
    tree->SetBranchAddress("ED_mcs_seg_fit_Z2", _mcs_seg_fit_Z2);
  } else {
    _nMcsSegments = 0;
    _mcs_mainTrack_uniqueID = -1;
  }

  _edHasMcsTrueBranches = tree->GetBranch("ED_nMcsTrueSegments") != nullptr;
  if (_edHasMcsTrueBranches) {
    tree->SetBranchAddress("ED_nMcsTrueSegments", &_nMcsTrueSegments);
    tree->SetBranchAddress("ED_mcs_true_traj_trueID", &_mcs_true_traj_trueID);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_X1", _mcs_true_seg_chord_X1);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_Y1", _mcs_true_seg_chord_Y1);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_Z1", _mcs_true_seg_chord_Z1);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_X2", _mcs_true_seg_chord_X2);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_Y2", _mcs_true_seg_chord_Y2);
    tree->SetBranchAddress("ED_mcs_true_seg_chord_Z2", _mcs_true_seg_chord_Z2);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_X1", _mcs_true_seg_fit_X1);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_Y1", _mcs_true_seg_fit_Y1);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_Z1", _mcs_true_seg_fit_Z1);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_X2", _mcs_true_seg_fit_X2);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_Y2", _mcs_true_seg_fit_Y2);
    tree->SetBranchAddress("ED_mcs_true_seg_fit_Z2", _mcs_true_seg_fit_Z2);
  } else {
    _nMcsTrueSegments = 0;
    _mcs_true_traj_trueID = -999;
  }
  return true;
}

//********************************************************************
void pionMomentumEventDisplay::DrawParticles3D(TEveScene* scene) {
//********************************************************************
  pdEventDisplay::DrawParticles3D(scene);
  if (!scene) return;

  TEveElement* particlesElement = nullptr;
  for (TEveElement::List_i it = scene->BeginChildren(); it != scene->EndChildren(); ++it) {
    TEveElement* child = *it;
    if (!child) continue;
    const char* name = child->GetElementName();
    if (name && std::string(name) == "Particles") {
      // Keep the last "Particles" group, which corresponds to the latest draw call.
      particlesElement = child;
    }
  }
  if (!particlesElement) return;

  Int_t trjIndex = 0;
  Int_t particleSlot = 0;
  for (TEveElement::List_i it = particlesElement->BeginChildren(); it != particlesElement->EndChildren(); ++it) {
    if (particleSlot >= _nParticles || particleSlot >= kMaxParticles) break;

    TEveElementList* particleGroup = dynamic_cast<TEveElementList*>(*it);
    if (!particleGroup) {
      ++particleSlot;
      continue;
    }

    const Int_t nTrj = std::max(0, _particle_nTrjPoints[particleSlot]);
    TEveElementList* trjGroup = new TEveElementList("Reco Trajectory Points (reco-SCE)");
    particleGroup->AddElement(trjGroup);

    if (nTrj > 0) {
      Int_t color = GetParticleColor(_particle_PDG[particleSlot]);
      if (color == kBlack) color = kGray + 1;

      TEvePointSet* firstTrjSet = new TEvePointSet("First Reco TrjPoint (reco-SCE)");
      firstTrjSet->SetMarkerStyle(29);
      firstTrjSet->SetMarkerSize(2.0);
      firstTrjSet->SetMainColor(color);
      trjGroup->AddElement(firstTrjSet);

      TEvePointSet* trjSet = new TEvePointSet("Reco TrjPoints (reco-SCE)");
      trjSet->SetMarkerStyle(20);
      trjSet->SetMarkerSize(0.55);
      trjSet->SetMainColor(color);
      trjGroup->AddElement(trjSet);

      for (Int_t t = 0; t < nTrj && trjIndex < _nTrjPoints && trjIndex < kMaxTrjPoints; ++t) {
        const Float_t x = _trj_X[trjIndex];
        const Float_t y = _trj_Y[trjIndex];
        const Float_t z = _trj_Z[trjIndex];
        if (x > -900.f && y > -900.f && z > -900.f) {
          if (t == nTrj - 1) {
            firstTrjSet->SetNextPoint(x, y, z);
          } else {
            trjSet->SetNextPoint(x, y, z);
          }
        }
        ++trjIndex;
      }
    } else {
      trjIndex += nTrj;
    }

    ++particleSlot;
  }

  TEveElementList* mcsGroup = new TEveElementList("MCS");
  scene->AddElement(mcsGroup);

  if (_edHasMcsBranches && _nMcsSegments > 0) {
    TEveElementList* perReco = new TEveElementList(
        Form("Reco MainTrack UID=%d", _mcs_mainTrack_uniqueID));
    mcsGroup->AddElement(perReco);

    TEveElementList* chordGroup = new TEveElementList("MCS Segments (chord, reco-SCE)");
    TEveElementList* fitGroup   = new TEveElementList("MCS Segments (PCA fit, reco-SCE)");
    perReco->AddElement(chordGroup);
    perReco->AddElement(fitGroup);

    TEvePointSet* chordSegStarts = new TEvePointSet("Reco chord segment starts");
    TEvePointSet* chordSegEnds = new TEvePointSet("Reco chord segment ends");
    chordSegStarts->SetMarkerStyle(20);
    chordSegStarts->SetMarkerSize(0.9);
    chordSegStarts->SetMainColor(kRed);
    chordSegEnds->SetMarkerStyle(21);
    chordSegEnds->SetMarkerSize(0.9);
    chordSegEnds->SetMainColor(kBlue);

    TEvePointSet* fitSegStarts = new TEvePointSet("Reco PCA fit segment starts");
    TEvePointSet* fitSegEnds = new TEvePointSet("Reco PCA fit segment ends");
    fitSegStarts->SetMarkerStyle(24);
    fitSegStarts->SetMarkerSize(0.95);
    fitSegStarts->SetMainColor(kGreen + 2);
    fitSegEnds->SetMarkerStyle(25);
    fitSegEnds->SetMarkerSize(0.95);
    fitSegEnds->SetMainColor(kMagenta);

    const Int_t nSeg = std::min<Int_t>(_nMcsSegments, kMaxMcsSegments);
    for (Int_t s = 0; s < nSeg; ++s) {
      const Float_t cx1 = _mcs_seg_chord_X1[s];
      const Float_t cy1 = _mcs_seg_chord_Y1[s];
      const Float_t cz1 = _mcs_seg_chord_Z1[s];
      const Float_t cx2 = _mcs_seg_chord_X2[s];
      const Float_t cy2 = _mcs_seg_chord_Y2[s];
      const Float_t cz2 = _mcs_seg_chord_Z2[s];
      if (cx1 > -900.f && cy1 > -900.f && cz1 > -900.f &&
          cx2 > -900.f && cy2 > -900.f && cz2 > -900.f) {
        TEveLine* lc = new TEveLine(Form("reco_seg_chord_%d", s));
        lc->SetMainColor(kBlack);
        lc->SetLineWidth(2);
        lc->SetLineStyle(1);
        lc->SetNextPoint(cx1, cy1, cz1);
        lc->SetNextPoint(cx2, cy2, cz2);
        chordGroup->AddElement(lc);
        chordSegStarts->SetNextPoint(cx1, cy1, cz1);
        chordSegEnds->SetNextPoint(cx2, cy2, cz2);
      }

      const Float_t fx1 = _mcs_seg_fit_X1[s];
      const Float_t fy1 = _mcs_seg_fit_Y1[s];
      const Float_t fz1 = _mcs_seg_fit_Z1[s];
      const Float_t fx2 = _mcs_seg_fit_X2[s];
      const Float_t fy2 = _mcs_seg_fit_Y2[s];
      const Float_t fz2 = _mcs_seg_fit_Z2[s];
      if (fx1 > -900.f && fy1 > -900.f && fz1 > -900.f &&
          fx2 > -900.f && fy2 > -900.f && fz2 > -900.f) {
        TEveLine* lf = new TEveLine(Form("reco_seg_fit_%d", s));
        lf->SetMainColor(kBlack);
        lf->SetLineWidth(2);
        lf->SetLineStyle(1);
        lf->SetNextPoint(fx1, fy1, fz1);
        lf->SetNextPoint(fx2, fy2, fz2);
        fitGroup->AddElement(lf);
        fitSegStarts->SetNextPoint(fx1, fy1, fz1);
        fitSegEnds->SetNextPoint(fx2, fy2, fz2);
      }
    }

    chordGroup->AddElement(chordSegStarts);
    chordGroup->AddElement(chordSegEnds);
    fitGroup->AddElement(fitSegStarts);
    fitGroup->AddElement(fitSegEnds);
  }

  if (_edHasMcsTrueBranches && _nMcsTrueSegments > 0) {
    TEveElementList* perTrue = new TEveElementList(
        Form("True trajectory TrueID=%d (SCE-distorted G4)", _mcs_true_traj_trueID));
    mcsGroup->AddElement(perTrue);

    TEveElementList* trueChordGroup = new TEveElementList("MCS Segments (chord, true SCE-distorted G4)");
    TEveElementList* trueFitGroup   = new TEveElementList("MCS Segments (PCA fit, true SCE-distorted G4)");
    perTrue->AddElement(trueChordGroup);
    perTrue->AddElement(trueFitGroup);

    TEvePointSet* trueChordStarts = new TEvePointSet("True chord segment starts");
    TEvePointSet* trueChordEnds = new TEvePointSet("True chord segment ends");
    trueChordStarts->SetMarkerStyle(20);
    trueChordStarts->SetMarkerSize(0.9);
    trueChordStarts->SetMainColor(kRed);
    trueChordEnds->SetMarkerStyle(21);
    trueChordEnds->SetMarkerSize(0.9);
    trueChordEnds->SetMainColor(kBlue);

    TEvePointSet* trueFitStarts = new TEvePointSet("True PCA fit segment starts");
    TEvePointSet* trueFitEnds = new TEvePointSet("True PCA fit segment ends");
    trueFitStarts->SetMarkerStyle(24);
    trueFitStarts->SetMarkerSize(0.95);
    trueFitStarts->SetMainColor(kGreen + 2);
    trueFitEnds->SetMarkerStyle(25);
    trueFitEnds->SetMarkerSize(0.95);
    trueFitEnds->SetMainColor(kMagenta);

    const Int_t nTrueSeg = std::min<Int_t>(_nMcsTrueSegments, kMaxMcsSegments);
    for (Int_t s = 0; s < nTrueSeg; ++s) {
      const Float_t cx1 = _mcs_true_seg_chord_X1[s];
      const Float_t cy1 = _mcs_true_seg_chord_Y1[s];
      const Float_t cz1 = _mcs_true_seg_chord_Z1[s];
      const Float_t cx2 = _mcs_true_seg_chord_X2[s];
      const Float_t cy2 = _mcs_true_seg_chord_Y2[s];
      const Float_t cz2 = _mcs_true_seg_chord_Z2[s];
      if (cx1 > -900.f && cy1 > -900.f && cz1 > -900.f &&
          cx2 > -900.f && cy2 > -900.f && cz2 > -900.f) {
        TEveLine* lc = new TEveLine(Form("true_seg_chord_%d", s));
        lc->SetMainColor(kBlack);
        lc->SetLineWidth(2);
        lc->SetLineStyle(1);
        lc->SetNextPoint(cx1, cy1, cz1);
        lc->SetNextPoint(cx2, cy2, cz2);
        trueChordGroup->AddElement(lc);
        trueChordStarts->SetNextPoint(cx1, cy1, cz1);
        trueChordEnds->SetNextPoint(cx2, cy2, cz2);
      }

      const Float_t fx1 = _mcs_true_seg_fit_X1[s];
      const Float_t fy1 = _mcs_true_seg_fit_Y1[s];
      const Float_t fz1 = _mcs_true_seg_fit_Z1[s];
      const Float_t fx2 = _mcs_true_seg_fit_X2[s];
      const Float_t fy2 = _mcs_true_seg_fit_Y2[s];
      const Float_t fz2 = _mcs_true_seg_fit_Z2[s];
      if (fx1 > -900.f && fy1 > -900.f && fz1 > -900.f &&
          fx2 > -900.f && fy2 > -900.f && fz2 > -900.f) {
        TEveLine* lf = new TEveLine(Form("true_seg_fit_%d", s));
        lf->SetMainColor(kBlack);
        lf->SetLineWidth(2);
        lf->SetLineStyle(1);
        lf->SetNextPoint(fx1, fy1, fz1);
        lf->SetNextPoint(fx2, fy2, fz2);
        trueFitGroup->AddElement(lf);
        trueFitStarts->SetNextPoint(fx1, fy1, fz1);
        trueFitEnds->SetNextPoint(fx2, fy2, fz2);
      }
    }

    trueChordGroup->AddElement(trueChordStarts);
    trueChordGroup->AddElement(trueChordEnds);
    trueFitGroup->AddElement(trueFitStarts);
    trueFitGroup->AddElement(trueFitEnds);
  }
}
