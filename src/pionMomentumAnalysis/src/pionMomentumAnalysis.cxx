#include "pionMomentumAnalysis.hxx"
#include "pionMomentumSelection.hxx"
#include "pionMomentumTree.hxx"
#include "pionMomentumAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"
#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"
#include "pdMomReconstruction.hxx"
#include "pdUtilsDEdx.hxx"
#include "pdUtilsRangePID.hxx"
#include "pdUtilsTrack.hxx"
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

//********************************************************************
pionMomentumAnalysis::pionMomentumAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************
  _ApplySCECorrection = false;
  _ApplySCESystematic = false;
  _MCSRadiationLengthCm = 14.0;
  _MCSTargetSegmentLengthCm = 10.0;
  _MCSMinSegmentLengthCm = 0.5;
  _MCStheta0FloorRad = 1e-6;
  _MCSMaxAbsDeltaThetaRad = -1.0;
  _MCSDropFirstNSegments = 0;
  _MCSDropLastNSegments = 0;
  _MCSUseDetectorSigma = false;
  _MCSDetectorSigmaFloorRad = 1e-6;
  _MCSDetectorSigmaA = 0.0;
  _MCSDetectorSigmaC = 0.0;
  _MCSMomentumScanMaxGeV = 2.5;
  _TLEMinInteriorHits = 15;
  _TLESkipHitsFirst = 3;
  _TLESkipHitsLast = 3;
  _TLEDedxMinMeVcm = 0.5;
  _TLEDedxMaxMeVcm = 5.0;
  _TLEScanLmaxCm = 450.0;
  _TLEScanStepCm = 1.0;
  _TLEScanStepFineCm = 0.25;
  _TLELowPMomentumRefineGeV = 0.4;
  _TLEDedxPdfPathCm = 0.65;
  _BraggChi2PiMaxResidualRangeCm = 26.0;
  _BraggChi2PiSigmaMeVcm = 0.5;
  _BraggChi2PiMinHits = 3;
  _BraggChi2PiSkipHitsFirst = 0;
  _BraggChi2PiSkipHitsLast = 0;
  _BraggChi2PiDedxMinMeVcm = 0.;
  _BraggChi2PiDedxMaxMeVcm = 0.;
  _MomentumRequireTrueBeamPion = false;
  _RunTLEFit = true;
  _RunMCS = true;
  _RunCSDA = true;
  _RunBeamDaughterTLETruncationScan = false;
  _RunBeamDaughterMCSTruncationScan = false;
  _EnableMomentumDiagnosticMultigraphs = false;
  _FreeRangeComputeDedxBiasDiagnostics = false;
  _StoppingMaxTrueEndMomentumGeV = 0.03;
  _CreateEventDisplay = false;
  _EventDisplay = nullptr;
}

//********************************************************************
pionMomentumAnalysis::~pionMomentumAnalysis() {
//********************************************************************
  if (_EventDisplay) {
    delete _EventDisplay;
    _EventDisplay = nullptr;
  }
}

//********************************************************************
bool pionMomentumAnalysis::Initialize() {
//********************************************************************
  if (!baseAnalysis::Initialize()) return false;

  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pionMomentumAnalysis.MinAccumLevelToSave"));
  _ApplySCECorrection = ND::params().GetParameterI("pionMomentumAnalysis.ApplySCECorrection");
  _ApplySCESystematic = ND::params().GetParameterI("pionMomentumAnalysis.ApplySCESystematic");
  {
    pdMomReconstruction::PionMCSConfig mcsLoaded;
    pdMomReconstruction::FillPionMCSConfig_FromPionMomentumParams(mcsLoaded);
    _MCSRadiationLengthCm = mcsLoaded.likelihood.radiationLengthCm;
    _MCSTargetSegmentLengthCm = mcsLoaded.likelihood.targetSegmentLengthCm;
    _MCSMinSegmentLengthCm = mcsLoaded.likelihood.minSegmentLengthCm;
    _MCStheta0FloorRad = mcsLoaded.likelihood.theta0FloorRad;
    _MCSMaxAbsDeltaThetaRad = mcsLoaded.likelihood.maxAbsDeltaThetaRad;
    _MCSDropFirstNSegments = mcsLoaded.dropFirstNSegments;
    _MCSDropLastNSegments = mcsLoaded.dropLastNSegments;
    _MCSUseDetectorSigma = mcsLoaded.likelihood.useDetectorSigma;
    _MCSDetectorSigmaFloorRad = mcsLoaded.likelihood.detectorSigmaFloorRad;
    _MCSDetectorSigmaA = mcsLoaded.likelihood.detectorSigmaA;
    _MCSDetectorSigmaC = mcsLoaded.likelihood.detectorSigmaC;
    _MCSMomentumScanMaxGeV = mcsLoaded.momentumScanMaxGeV;
  }
  {
    pdMomReconstruction::PionTLEFitConfig tleLoaded;
    pdMomReconstruction::FillPionTLEFitConfig_FromPionMomentumParams(tleLoaded);
    _TLEMinInteriorHits = tleLoaded.minInteriorHits;
    _TLESkipHitsFirst = tleLoaded.skipHitsFirst;
    _TLESkipHitsLast = tleLoaded.skipHitsLast;
    _TLEDedxMinMeVcm = tleLoaded.dedxMinMeVcm;
    _TLEDedxMaxMeVcm = tleLoaded.dedxMaxMeVcm;
    _TLEScanLmaxCm = tleLoaded.scanLmaxCm;
    _TLEScanStepCm = tleLoaded.scanStepCm;
    _TLEScanStepFineCm = tleLoaded.scanStepFineCm;
    _TLELowPMomentumRefineGeV = tleLoaded.lowPMomentumRefineGeV;
    _TLEDedxPdfPathCm = tleLoaded.dedxPdfPathCm;
  }
  _BraggChi2PiMaxResidualRangeCm = ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiMaxResidualRangeCm");
  _BraggChi2PiSigmaMeVcm = ND::params().HasParameter("pionMomentumAnalysis.BraggChi2PiSigmaMeVcm")
                               ? ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiSigmaMeVcm")
                               : 0.5;
  _BraggChi2PiMinHits = ND::params().GetParameterI("pionMomentumAnalysis.BraggChi2PiMinHits");
  _BraggChi2PiSkipHitsFirst = ND::params().GetParameterI("pionMomentumAnalysis.BraggChi2PiSkipHitsFirst");
  _BraggChi2PiSkipHitsLast = ND::params().GetParameterI("pionMomentumAnalysis.BraggChi2PiSkipHitsLast");
  _BraggChi2PiDedxMinMeVcm = ND::params().HasParameter("pionMomentumAnalysis.BraggChi2PiDedxMinMeVcm")
                                ? ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiDedxMinMeVcm")
                                : 0.;
  _BraggChi2PiDedxMaxMeVcm = ND::params().HasParameter("pionMomentumAnalysis.BraggChi2PiDedxMaxMeVcm")
                                ? ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiDedxMaxMeVcm")
                                : 0.;
  _MomentumRequireTrueBeamPion = ND::params().GetParameterI("pionMomentumAnalysis.MomentumRequireTrueBeamPion");
  _RunTLEFit = ND::params().GetParameterI("pionMomentumAnalysis.RunTLEFit");
  _RunMCS = ND::params().GetParameterI("pionMomentumAnalysis.RunMCS");
  _RunCSDA = ND::params().GetParameterI("pionMomentumAnalysis.RunCSDA");
  _RunBeamDaughterTLETruncationScan =
      ND::params().HasParameter("pionMomentumAnalysis.RunBeamDaughterTLETruncationScan")
          ? ND::params().GetParameterI("pionMomentumAnalysis.RunBeamDaughterTLETruncationScan")
          : 0;
  _RunBeamDaughterMCSTruncationScan =
      ND::params().HasParameter("pionMomentumAnalysis.RunBeamDaughterMCSTruncationScan")
          ? ND::params().GetParameterI("pionMomentumAnalysis.RunBeamDaughterMCSTruncationScan")
          : 0;
  _EnableMomentumDiagnosticMultigraphs =
      ND::params().GetParameterI("pionMomentumAnalysis.EnableMomentumDiagnosticMultigraphs");
  _FreeRangeComputeDedxBiasDiagnostics =
      ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeComputeDedxBiasDiagnostics");
  _StoppingMaxTrueEndMomentumGeV = ND::params().GetParameterD("pionMomentumAnalysis.StoppingMaxTrueEndMomentumGeV");
  _CreateEventDisplay = ND::params().GetParameterI("pionMomentumAnalysis.CreateEventDisplay");

  // Register the standard category set used by anaUtils::FillCategories(mainTrack, "")
  anaUtils::AddStandardCategories();
  pionMomentumAnaUtils::AddCustomCategories();

  if (_CreateEventDisplay && !_EventDisplay) {
    _EventDisplay = new pionMomentumEventDisplay();
  }

  return true;
}

//********************************************************************
void pionMomentumAnalysis::Finalize() {
//********************************************************************
  pionMomentumAnaUtils::WriteMomentumDiagnostics(output());
}

//********************************************************************
void pionMomentumAnalysis::DefineInputConverters() {
//********************************************************************
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void pionMomentumAnalysis::DefineSelections() {
//********************************************************************
  sel().AddSelection("pionMomentumSelection", "Pion momentum analysis (pass all)", new pionMomentumSelection(false));
}

//********************************************************************
void pionMomentumAnalysis::DefineCorrections() {
//********************************************************************
  baseAnalysis::DefineCorrections();
  if (_ApplySCECorrection)
    corr().AddCorrection(0, "sce geometric correction", new ParticlePositionSCECorrection());
}

//********************************************************************
void pionMomentumAnalysis::DefineSystematics() {
//********************************************************************
  baseAnalysis::DefineSystematics();
  if (_ApplySCESystematic)
    evar().AddEventVariation(kSCEGeometric, "SCE variation", new SCEGeometricVariation());
}

//********************************************************************
void pionMomentumAnalysis::DefineConfigurations() {
//********************************************************************
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void pionMomentumAnalysis::DefineMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
  pionMomentumTree::AddPionMomentumVariables_BeamParticleReco(output());
  pionMomentumTree::AddPionMomentumVariables_BeamTrueDaughterBragg(output());
  pionMomentumTree::AddPionMomentumVariables_BeamDaughterTleTruncScan(output());
  pionMomentumTree::AddPionMomentumVariables_BeamDaughterMcsTruncScan(output());

  if (_CreateEventDisplay && _EventDisplay) {
    _EventDisplay->InitializeTree(output());
  }
}

//********************************************************************
void pionMomentumAnalysis::DefineTruthTree() {
//********************************************************************
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void pionMomentumAnalysis::FillMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  std::vector<double> mainthetascatter;
  std::vector<double> mainthetascatterxz;
  std::vector<double> mainthetascatteryz;
  std::vector<double> mainsegmentlengthpair;
  std::vector<double> mainsegmentlengthsingle;
  std::vector<double> mainsegmentfitchi2ndfsingle;
  std::vector<double> mainsegmentrrtoendsingle;
  std::vector<double> mainthetascatterused;
  std::vector<double> mainthetascatterxzusd;
  std::vector<double> mainthetascatteryzusd;
  std::vector<double> mainrrscatterused;
  double mainmcsmomentum = -999.0;
  double maintlefitmomentum = -999.0;
  std::vector<double> pAxisTleGeV;
  std::vector<double> logLTle;
  std::vector<double> pAxisMcsGeV;
  std::vector<double> logLMcs;
  std::vector<double> maintruethetascatter;
  std::vector<double> maintruethetascatterxz;
  std::vector<double> maintruethetascatteryz;
  std::vector<double> maintruesegmentlengthpair;
  std::vector<double> maintruesegmentlengthsingle;
  std::vector<double> maintruesegmentfitchi2ndfsingle;
  std::vector<double> maintruesegmentrrtoendsingle;
  double maintruemcsmomentum = -999.0;
  Int_t mainntruetrajpointsTree = -999;
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
  AnaParticlePD* beampart = NULL;
  if (beam && beam->BeamParticle) beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);

  if (box().MainTrack) {
    bool trueBeamPion = false;
    if (box().MainTrack->TrueObject && beampart && beampart->TrueObject) {
      const AnaTrueParticlePD* mainTrue = static_cast<const AnaTrueParticlePD*>(box().MainTrack->TrueObject);
      const AnaTrueParticlePD* beamTrue = static_cast<const AnaTrueParticlePD*>(beampart->TrueObject);
      if (mainTrue && beamTrue) {
        trueBeamPion = (std::abs(mainTrue->PDG) == 211) && (mainTrue->ID == beamTrue->ID);
      }
    }
    const bool runMomentum = !_MomentumRequireTrueBeamPion || trueBeamPion;
    const bool runTLE = runMomentum && _RunTLEFit;
    const bool runMCS = runMomentum && _RunMCS;

    // True trajectory for MCS: beam TrueObject gets TrjPoints in PDSP converter; MainTrack->TrueObject can be same ID
    // but a different clone with empty TrjPoints — prefer beam when IDs match and main has no points.
    const AnaTrueParticlePD* mainTrueParticle =
        box().MainTrack->TrueObject ? static_cast<const AnaTrueParticlePD*>(box().MainTrack->TrueObject) : nullptr;
    const AnaTrueParticlePD* beamTrueParticle =
        (beampart && beampart->TrueObject) ? static_cast<const AnaTrueParticlePD*>(beampart->TrueObject) : nullptr;
    const AnaTrueParticlePD* trajTruthForTrueMcs = nullptr;
    if (mainTrueParticle && beamTrueParticle && mainTrueParticle->ID == beamTrueParticle->ID &&
        mainTrueParticle->TrjPoints.empty() && !beamTrueParticle->TrjPoints.empty()) {
      trajTruthForTrueMcs = beamTrueParticle;
    } else if (mainTrueParticle) {
      trajTruthForTrueMcs = mainTrueParticle;
    } else {
      trajTruthForTrueMcs = beamTrueParticle;
    }
    // Beam PFP truth (byHits/byE) can be a different ID than true_beam_ID with no step points;
    // instrument beam truth still gets true_beam_traj_* in PDSPAnalyzerTreeConverter.
    if (trajTruthForTrueMcs && trajTruthForTrueMcs->TrjPoints.empty() && beamTrueParticle &&
        !beamTrueParticle->TrjPoints.empty()) {
      trajTruthForTrueMcs = beamTrueParticle;
    }
    if (trajTruthForTrueMcs) mainntruetrajpointsTree = static_cast<Int_t>(trajTruthForTrueMcs->TrjPoints.size());

    if (const char* dbgTrj = std::getenv("HIGHLAND_DEBUG_TRUE_TRJ_COUNTS");
        dbgTrj && std::strcmp(dbgTrj, "1") == 0 && GetEvent().EventInfo) {
      const AnaEventInfoB& ei = *GetEvent().EventInfo;
      const int mainId = mainTrueParticle ? mainTrueParticle->ID : -999;
      const int beamId = beamTrueParticle ? beamTrueParticle->ID : -999;
      const int trajId = trajTruthForTrueMcs ? trajTruthForTrueMcs->ID : -999;
      const size_t nMain = mainTrueParticle ? mainTrueParticle->TrjPoints.size() : 0;
      const size_t nBeam = beamTrueParticle ? beamTrueParticle->TrjPoints.size() : 0;
      std::cerr << "[TRUE_TRJ_COUNTS] stage=analysis_ntuple_truth run=" << ei.Run << " subrun=" << ei.SubRun
                << " evt=" << ei.Event << " mainTrueID=" << mainId << " beamTrueID=" << beamId
                << " trajUsedID=" << trajId << " nTrj_mainTrueObject=" << nMain
                << " nTrj_beamTrueObject=" << nBeam << " nTrj_usedForTree=" << mainntruetrajpointsTree << "\n";
    }

    if (runTLE) {
      pdMomReconstruction::PionTLEFitConfig tleCfg;
      pdMomReconstruction::FillPionTLEFitConfig_FromPionMomentumParams(tleCfg);
      pdMomReconstruction::PionTLEFitResult tleOut;
      const bool wantTleCurve = _EnableMomentumDiagnosticMultigraphs;
      if (pdMomReconstruction::EstimatePionMomentumTLEFit(box().MainTrack, tleCfg, tleOut,
                                                        wantTleCurve ? &pAxisTleGeV : nullptr,
                                                        wantTleCurve ? &logLTle : nullptr) &&
          tleOut.valid)
        maintlefitmomentum = tleOut.bestMomentumGeV;
    }

    if (runMCS) {
      pdMomReconstruction::PionMCSConfig mcsCfg;
      pdMomReconstruction::FillPionMCSConfig_FromPionMomentumParams(mcsCfg);

      pdMomReconstruction::PionMCSRecoBuffers recoBufs;
      pdMomReconstruction::PionMCSResult recoMom;
      const bool wantMcsCurve = _EnableMomentumDiagnosticMultigraphs;
      const bool recoSamplesOk = pdMomReconstruction::EstimatePionMomentumMCSReco(
          box().MainTrack, mcsCfg, recoMom, recoBufs, wantMcsCurve ? &pAxisMcsGeV : nullptr,
          wantMcsCurve ? &logLMcs : nullptr);
      if (recoSamplesOk) {
        mainthetascatter = std::move(recoBufs.deltaThetaPair);
        mainthetascatterxz = std::move(recoBufs.deltaThetaXz);
        mainthetascatteryz = std::move(recoBufs.deltaThetaYz);
        mainsegmentlengthpair = std::move(recoBufs.segmentLengthPair);
        mainsegmentlengthsingle = std::move(recoBufs.segmentLengthSingle);
        mainsegmentfitchi2ndfsingle = std::move(recoBufs.segmentFitChi2NdfSingle);
        mainsegmentrrtoendsingle = std::move(recoBufs.segmentRrToEndSingle);
        mainthetascatterused = std::move(recoBufs.scatterUsedPair);
        mainthetascatterxzusd = std::move(recoBufs.scatterUsedXz);
        mainthetascatteryzusd = std::move(recoBufs.scatterUsedYz);
        mainrrscatterused = std::move(recoBufs.scatterUsedRR);
        mainmcsmomentum = recoMom.valid ? recoMom.bestMomentumGeV : -999.0;
      } else {
        mainthetascatter.clear();
        mainthetascatterxz.clear();
        mainthetascatteryz.clear();
        mainsegmentlengthpair.clear();
        mainsegmentlengthsingle.clear();
        mainsegmentfitchi2ndfsingle.clear();
        mainsegmentrrtoendsingle.clear();
        mainthetascatterused.clear();
        mainthetascatterxzusd.clear();
        mainthetascatteryzusd.clear();
        mainrrscatterused.clear();
        mainmcsmomentum = -999.0;
      }

      // True-trajectory MCS: same cfg and drop logic as reco; positions = AnaTrueTrajectoryPointPD::Position (SCE).
      std::vector<TVector3> trueOrdered;
      if (trajTruthForTrueMcs) {
        trueOrdered.reserve(trajTruthForTrueMcs->TrjPoints.size());
        for (const auto& tp : trajTruthForTrueMcs->TrjPoints) {
          if (!tp.IsInTPC) continue;
          const TVector3& p = tp.Position;
          if (std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900.0 &&
              p.Y() > -900.0 && p.Z() > -900.0)
            trueOrdered.push_back(p);
        }
      }

      pdMomReconstruction::PionMCSTrueBuffers trueBufs;
      pdMomReconstruction::PionMCSResult trueMom;
      const bool trueSamplesOk =
          pdMomReconstruction::EstimatePionMomentumMCSTrue(trueOrdered, mcsCfg, trueMom, trueBufs);
      if (trueSamplesOk) {
        maintruethetascatter = std::move(trueBufs.deltaThetaPair);
        maintruethetascatterxz = std::move(trueBufs.deltaThetaXz);
        maintruethetascatteryz = std::move(trueBufs.deltaThetaYz);
        maintruesegmentlengthpair = std::move(trueBufs.segmentLengthPair);
        maintruesegmentlengthsingle = std::move(trueBufs.segmentLengthSingle);
        maintruesegmentfitchi2ndfsingle = std::move(trueBufs.segmentFitChi2NdfSingle);
        maintruesegmentrrtoendsingle = std::move(trueBufs.segmentRrToEndSingle);
        maintruemcsmomentum = trueMom.valid ? trueMom.bestMomentumGeV : -999.0;
      } else {
        maintruethetascatter.clear();
        maintruethetascatterxz.clear();
        maintruethetascatteryz.clear();
        maintruesegmentlengthpair.clear();
        maintruesegmentlengthsingle.clear();
        maintruesegmentfitchi2ndfsingle.clear();
        maintruesegmentrrtoendsingle.clear();
        maintruemcsmomentum = -999.0;
      }

    } else {
      mainthetascatter.clear();
      mainthetascatterxz.clear();
      mainthetascatteryz.clear();
      mainsegmentlengthpair.clear();
      mainsegmentlengthsingle.clear();
      mainsegmentfitchi2ndfsingle.clear();
      mainsegmentrrtoendsingle.clear();
      mainthetascatterused.clear();
      mainthetascatterxzusd.clear();
      mainthetascatteryzusd.clear();
      mainrrscatterused.clear();
      mainmcsmomentum = -999.0;
      maintruethetascatter.clear();
      maintruethetascatterxz.clear();
      maintruethetascatteryz.clear();
      maintruesegmentlengthpair.clear();
      maintruesegmentlengthsingle.clear();
      maintruesegmentfitchi2ndfsingle.clear();
      maintruesegmentrrtoendsingle.clear();
      maintruemcsmomentum = -999.0;
      pAxisMcsGeV.clear();
      logLMcs.clear();
    }

    const int stoppingCodeForDiag =
        pionMomentumAnaUtils::MainTrueStoppingCode(box().MainTrack, _StoppingMaxTrueEndMomentumGeV);
    if (_EnableMomentumDiagnosticMultigraphs && runMomentum) {
    pionMomentumAnaUtils::MomentumDiagConfig diagCfg;
    diagCfg.enableMomentumDiagnosticMultigraphs = _EnableMomentumDiagnosticMultigraphs;
    diagCfg.freeRangeComputeDedxBiasDiagnostics = _FreeRangeComputeDedxBiasDiagnostics;
    diagCfg.freeRangeDedxMinInteriorHits = _TLEMinInteriorHits;
    diagCfg.freeRangeDedxSkipHitsFirst = _TLESkipHitsFirst;
    diagCfg.freeRangeDedxSkipHitsLast = _TLESkipHitsLast;
    diagCfg.freeRangeDedxMinMeVcm = _TLEDedxMinMeVcm;
    diagCfg.freeRangeDedxMaxMeVcm = _TLEDedxMaxMeVcm;
    diagCfg.freeRangeScanLmaxCm = _TLEScanLmaxCm;
    diagCfg.freeRangeScanStepCm = _TLEScanStepCm;
    diagCfg.freeRangeDedxPdfPathCm = _TLEDedxPdfPathCm;
    diagCfg.mcsRadiationLengthCm = _MCSRadiationLengthCm;
    const bool trueMcsDiagOk =
        !maintruethetascatter.empty() && maintruethetascatter.size() == maintruesegmentlengthpair.size();
    pionMomentumAnaUtils::MaybeAccumulateMainTrackMomentumDiagnostics(
        box().MainTrack, GetEvent(), stoppingCodeForDiag, runTLE, diagCfg, maintlefitmomentum, mainmcsmomentum, pAxisTleGeV,
        logLTle, pAxisMcsGeV, logLMcs, mainthetascatter, mainsegmentlengthpair, mainthetascatterused,
        mainrrscatterused, trueMcsDiagOk ? &maintruethetascatter : nullptr,
        trueMcsDiagOk ? &maintruesegmentlengthpair : nullptr);
    }
  }

  // Bragg-window mean χ²/hit: same construction as Chi2PID(211) (range–dE/dx template + errors), same RR/dE/dx/skip
  // window as former Eq. 6.1 branch. Runs whenever MainTrack exists — not gated on RunTLEFit or RunMCS.
  double mainbraggchi2pibb = -999.;
  Int_t mainbraggdedxnhits = -999;
  if (box().MainTrack) {
    double chi2 = -999.;
    int nh = 0;
    if (pdAnaUtils::ComputePionBraggWindowChi2PionRangeTemplate(
            box().MainTrack, _BraggChi2PiMaxResidualRangeCm, _BraggChi2PiMinHits, _BraggChi2PiSkipHitsFirst,
            _BraggChi2PiSkipHitsLast, _BraggChi2PiDedxMinMeVcm, _BraggChi2PiDedxMaxMeVcm, chi2, nh)) {
      mainbraggchi2pibb = chi2;
      mainbraggdedxnhits = static_cast<Int_t>(nh);
    }
  }

  std::vector<Int_t> beamDaughterTrueIds;
  std::vector<Int_t> beamDaughterTruePdgs;
  std::vector<Float_t> beamDaughterTrueP0GeV;
  std::vector<Float_t> beamDaughterTruePendGeV;
  std::vector<Float_t> beamDaughterTrueEkin0GeV;
  std::vector<Float_t> beamDaughterTrueEkinEndGeV;
  std::vector<Float_t> beamDaughterCSDAMomentumGeV;
  std::vector<Float_t> beamDaughterCSDAKineticEnergyGeV;
  std::vector<Int_t> beamDaughterTrueEndProcess;
  std::vector<Float_t> beamDaughterTrueEndToRecoStartDistCm;
  std::vector<Float_t> beamDaughterTrueEndToRecoEndDistCm;
  std::vector<double> beamDaughterBraggChi2;
  std::vector<Int_t> beamDaughterBraggNhits;
  std::vector<Int_t> beamDaughterTleTruncDauIdx;
  std::vector<Int_t> beamDaughterTleTruncK;
  std::vector<Int_t> beamDaughterTleTruncNhitsInt;
  std::vector<Float_t> beamDaughterTleTruncTrueEkin0GeV;
  std::vector<Float_t> beamDaughterTleTruncPtleGeV;
  std::vector<Float_t> beamDaughterTleTruncFracRes;
  std::vector<Int_t> beamDaughterMcsTruncDauIdx;
  std::vector<Int_t> beamDaughterMcsTruncK;
  std::vector<Int_t> beamDaughterMcsTruncNsegments;
  std::vector<Float_t> beamDaughterMcsTruncTrueEkin0GeV;
  std::vector<Float_t> beamDaughterMcsTruncPmcsGeV;
  std::vector<Float_t> beamDaughterMcsTruncFracRes;
  constexpr size_t kMaxBeamDaughtersTree = 64;
  constexpr size_t kMaxTleTruncRows = 4096;
  constexpr size_t kMaxMcsTruncRows = 4096;
  if (box().MainTrack) {
    bool trueBeamPionForDaughters = false;
    if (box().MainTrack->TrueObject && beampart && beampart->TrueObject) {
      const AnaTrueParticlePD* mainTrue = static_cast<const AnaTrueParticlePD*>(box().MainTrack->TrueObject);
      const AnaTrueParticlePD* beamTrue = static_cast<const AnaTrueParticlePD*>(beampart->TrueObject);
      if (mainTrue && beamTrue) {
        trueBeamPionForDaughters = (std::abs(mainTrue->PDG) == 211) && (mainTrue->ID == beamTrue->ID);
      }
    }
    const bool runMomentumForDaughters = !_MomentumRequireTrueBeamPion || trueBeamPionForDaughters;
    const bool runCSDA = runMomentumForDaughters && _RunCSDA;
    const std::vector<AnaRecObjectC*>& daughters = box().MainTrack->Daughters;
    for (const AnaRecObjectC* dauObj : daughters) {
      if (beamDaughterTrueIds.size() >= kMaxBeamDaughtersTree) break;
      const auto* recoDau = dynamic_cast<const AnaParticlePD*>(dauObj);
      if (!recoDau) continue;

      const Int_t outDauIdx = static_cast<Int_t>(beamDaughterTrueIds.size());

      Int_t trueId = -999;
      Int_t truePdg = -999;
      Float_t trueP0 = -999.f;
      Float_t truePend = -999.f;
      Float_t trueEkin0 = -999.f;
      Float_t trueEkinEnd = -999.f;
      Float_t csdaMomentum = -999.f;
      Float_t csdaKineticEnergy = -999.f;
      Int_t trueProcEnd = -999;
      Float_t trueEndToRecoStartCm = -999.f;
      Float_t trueEndToRecoEndCm = -999.f;
      if (recoDau->TrueObject) {
        const auto* t = static_cast<const AnaTrueParticlePD*>(recoDau->TrueObject);
        trueId = t ? t->ID : -999;
        truePdg = t ? t->PDG : -999;
        trueP0 = t ? static_cast<Float_t>(t->Momentum) : -999.f;
        truePend = t ? static_cast<Float_t>(t->MomentumEnd) : -999.f;
        if (t) {
          const double massMeV = pdMomReconstruction::GetRestMass(t->PDG);
          trueEkin0 =
              static_cast<Float_t>(pdMomReconstruction::MomentumToKineticEnergy(t->Momentum, massMeV) / 1000.0);
          trueEkinEnd =
              static_cast<Float_t>(pdMomReconstruction::MomentumToKineticEnergy(t->MomentumEnd, massMeV) / 1000.0);
        }
        trueProcEnd = t ? static_cast<Int_t>(t->ProcessEnd) : -999;
        if (t) {
          const TVector3 trueEnd(static_cast<double>(t->PositionEnd[0]), static_cast<double>(t->PositionEnd[1]),
                                 static_cast<double>(t->PositionEnd[2]));
          const TVector3 recoStart(static_cast<double>(recoDau->PositionStart[0]),
                                   static_cast<double>(recoDau->PositionStart[1]),
                                   static_cast<double>(recoDau->PositionStart[2]));
          const TVector3 recoEnd(static_cast<double>(recoDau->PositionEnd[0]),
                                 static_cast<double>(recoDau->PositionEnd[1]),
                                 static_cast<double>(recoDau->PositionEnd[2]));
          const auto validXYZ = [](const TVector3& p) {
            return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900.0 &&
                   p.Y() > -900.0 && p.Z() > -900.0;
          };
          if (validXYZ(trueEnd) && validXYZ(recoStart)) trueEndToRecoStartCm = static_cast<Float_t>((trueEnd - recoStart).Mag());
          if (validXYZ(trueEnd) && validXYZ(recoEnd)) trueEndToRecoEndCm = static_cast<Float_t>((trueEnd - recoEnd).Mag());
        }
      }
      beamDaughterTrueIds.push_back(trueId);
      beamDaughterTruePdgs.push_back(truePdg);
      beamDaughterTrueP0GeV.push_back(trueP0);
      beamDaughterTruePendGeV.push_back(truePend);
      beamDaughterTrueEkin0GeV.push_back(trueEkin0);
      beamDaughterTrueEkinEndGeV.push_back(trueEkinEnd);
      if (runCSDA) {
        csdaMomentum = pdMomReconstruction::EstimatePionMomentumFromCSDA(recoDau);
        csdaKineticEnergy = pdMomReconstruction::EstimatePionKineticEnergyFromCSDA(recoDau);
      }
      beamDaughterCSDAMomentumGeV.push_back(csdaMomentum);
      beamDaughterCSDAKineticEnergyGeV.push_back(csdaKineticEnergy);
      beamDaughterTrueEndProcess.push_back(trueProcEnd);
      beamDaughterTrueEndToRecoStartDistCm.push_back(trueEndToRecoStartCm);
      beamDaughterTrueEndToRecoEndDistCm.push_back(trueEndToRecoEndCm);

      double dChi2 = -999.;
      Int_t dNh = -999;
      if (recoDau->Type == AnaParticlePD::kTrack) {
        int nhTmp = 0;
        double chiTmp = -999.;
        AnaParticlePD* recoDauForBragg = const_cast<AnaParticlePD*>(recoDau);
        if (pdAnaUtils::ComputePionBraggWindowChi2PiEq61(
                recoDauForBragg, _BraggChi2PiMaxResidualRangeCm, _BraggChi2PiSigmaMeVcm,
                _BraggChi2PiMinHits, _BraggChi2PiSkipHitsFirst, _BraggChi2PiSkipHitsLast, _BraggChi2PiDedxMinMeVcm,
                _BraggChi2PiDedxMaxMeVcm, chiTmp, nhTmp)) {
          dChi2 = chiTmp;
          dNh = static_cast<Int_t>(nhTmp);
        }
      }
      beamDaughterBraggChi2.push_back(dChi2);
      beamDaughterBraggNhits.push_back(dNh);

      if (_RunBeamDaughterTLETruncationScan && runMomentumForDaughters && recoDau->Type == AnaParticlePD::kTrack) {
        const bool passPi = (std::abs(truePdg) == 211);
        const bool passChi2 = std::isfinite(dChi2) && dChi2 > -998. && dChi2 < 6.;
        const bool passEkin = (trueEkin0 > 0.f);
        if (passPi && passChi2 && passEkin) {
          pdMomReconstruction::PionTLEFitConfig tleCfg;
          pdMomReconstruction::FillPionTLEFitConfig_FromPionMomentumParams(tleCfg);
          AnaParticlePD* recoDauTle = const_cast<AnaParticlePD*>(recoDau);
          const double maxRR = std::numeric_limits<double>::max();
          const double noDedxWindowMin = 0.0;
          const double noDedxWindowMax = 0.0;
          std::vector<double> dedxProbe, rrProbe;
          if (pdAnaUtils::InteriorPionCollectionPlaneDedxRr(recoDauTle, maxRR, 0, tleCfg.skipHitsFirst,
                                                            tleCfg.skipHitsLast, noDedxWindowMin,
                                                            noDedxWindowMax, dedxProbe, rrProbe) &&
              dedxProbe.size() >= static_cast<size_t>(tleCfg.minInteriorHits)) {
            for (int k = 0;; ++k) {
              if (beamDaughterTleTruncDauIdx.size() >= kMaxTleTruncRows) break;
              std::vector<double> dedxI, rrI;
              if (!pdAnaUtils::InteriorPionCollectionPlaneDedxRr(recoDauTle, maxRR, 0, tleCfg.skipHitsFirst,
                                                                 tleCfg.skipHitsLast + k, noDedxWindowMin,
                                                                 noDedxWindowMax, dedxI, rrI))
                break;
              if (dedxI.size() < static_cast<size_t>(tleCfg.minInteriorHits)) break;
              pdMomReconstruction::PionTLEFitConfig cfgK = tleCfg;
              cfgK.skipHitsLast = tleCfg.skipHitsLast + k;
              pdMomReconstruction::PionTLEFitResult tleOut;
              if (pdMomReconstruction::EstimatePionMomentumTLEFit(recoDauTle, cfgK, tleOut) && tleOut.valid) {
                const double massPiMeV = pdMomReconstruction::GetRestMass(211);
                const double ekinTleMeV =
                    pdMomReconstruction::MomentumToKineticEnergy(tleOut.bestMomentumGeV, massPiMeV);
                const double ekinTrueMeV = static_cast<double>(trueEkin0) * 1000.0;
                if (ekinTrueMeV > 0.) {
                  const double fracRes = (ekinTleMeV - ekinTrueMeV) / ekinTrueMeV;
                  beamDaughterTleTruncDauIdx.push_back(outDauIdx);
                  beamDaughterTleTruncK.push_back(k);
                  beamDaughterTleTruncNhitsInt.push_back(static_cast<Int_t>(dedxI.size()));
                  beamDaughterTleTruncTrueEkin0GeV.push_back(trueEkin0);
                  beamDaughterTleTruncPtleGeV.push_back(static_cast<Float_t>(tleOut.bestMomentumGeV));
                  beamDaughterTleTruncFracRes.push_back(static_cast<Float_t>(fracRes));
                }
              }
            }
          }
        }
      }

      if (_RunBeamDaughterMCSTruncationScan && runMomentumForDaughters && recoDau->Type == AnaParticlePD::kTrack) {
        const bool passPi = (std::abs(truePdg) == 211);
        const bool passChi2 = std::isfinite(dChi2) && dChi2 > -998. && dChi2 < 6.;
        const bool passEkin = (trueEkin0 > 0.f);
        if (passPi && passChi2 && passEkin) {
          AnaParticlePD* recoDauMcs = const_cast<AnaParticlePD*>(recoDau);
          pdMomReconstruction::PionMCSConfig mcsCfgBase;
          pdMomReconstruction::FillPionMCSConfig_FromPionMomentumParams(mcsCfgBase);
          mcsCfgBase.dropFirstNSegments = 3;
          mcsCfgBase.dropLastNSegments = 3;
          for (int k = 0;; ++k) {
            if (beamDaughterMcsTruncDauIdx.size() >= kMaxMcsTruncRows) break;
            pdMomReconstruction::PionMCSConfig cfgK = mcsCfgBase;
            cfgK.dropLastNSegments = mcsCfgBase.dropLastNSegments + k;

            pdMomReconstruction::PionMCSResult mcsOut;
            pdMomReconstruction::PionMCSRecoBuffers mcsReco;
            if (!pdMomReconstruction::EstimatePionMomentumMCSReco(recoDauMcs, cfgK, mcsOut, mcsReco)) break;

            const int nMcs = static_cast<int>(mcsReco.deltaThetaPair.size());
            int dropFirst = std::max(0, cfgK.dropFirstNSegments);
            int dropLast = std::max(0, cfgK.dropLastNSegments);
            if (nMcs > 0) {
              dropFirst = std::min(dropFirst, nMcs - 1);
              dropLast = std::min(dropLast, nMcs - 1);
              if (dropFirst + dropLast >= nMcs) {
                dropFirst = 0;
                dropLast = 0;
              }
            }
            const int nSegmentsKept = std::max(0, nMcs - dropFirst - dropLast);
            if (nSegmentsKept < 15) break;

            if (mcsOut.valid) {
              const double massPiMeV = pdMomReconstruction::GetRestMass(211);
              const double ekinMcsMeV = pdMomReconstruction::MomentumToKineticEnergy(mcsOut.bestMomentumGeV, massPiMeV);
              const double ekinTrueMeV = static_cast<double>(trueEkin0) * 1000.0;
              if (ekinTrueMeV > 0.) {
                const double fracRes = (ekinMcsMeV - ekinTrueMeV) / ekinTrueMeV;
                beamDaughterMcsTruncDauIdx.push_back(outDauIdx);
                beamDaughterMcsTruncK.push_back(k);
                beamDaughterMcsTruncNsegments.push_back(nSegmentsKept);
                beamDaughterMcsTruncTrueEkin0GeV.push_back(trueEkin0);
                beamDaughterMcsTruncPmcsGeV.push_back(static_cast<Float_t>(mcsOut.bestMomentumGeV));
                beamDaughterMcsTruncFracRes.push_back(static_cast<Float_t>(fracRes));
              }
            }
          }
        }
      }
    }
  }

  pionMomentumTree::FillPionMomentumVariables_BeamParticleReco(
      output(), box().MainTrack, beampart, mainthetascatter, mainthetascatterxz, mainthetascatteryz,
      mainsegmentlengthsingle, mainsegmentfitchi2ndfsingle, mainsegmentrrtoendsingle, mainmcsmomentum,
      maintlefitmomentum, mainbraggchi2pibb, mainbraggdedxnhits, mainntruetrajpointsTree, maintruethetascatter,
      maintruethetascatterxz, maintruethetascatteryz, maintruesegmentlengthsingle, maintruesegmentfitchi2ndfsingle,
      maintruesegmentrrtoendsingle, maintruemcsmomentum);

  pionMomentumTree::FillPionMomentumVariables_BeamTrueDaughterBragg(
      output(), beamDaughterTrueIds, beamDaughterTruePdgs, beamDaughterTrueP0GeV, beamDaughterTruePendGeV,
      beamDaughterTrueEkin0GeV, beamDaughterTrueEkinEndGeV, beamDaughterCSDAMomentumGeV,
      beamDaughterCSDAKineticEnergyGeV, beamDaughterTrueEndProcess,
      beamDaughterTrueEndToRecoStartDistCm, beamDaughterTrueEndToRecoEndDistCm, beamDaughterBraggChi2,
      beamDaughterBraggNhits);

  pionMomentumTree::FillPionMomentumVariables_BeamDaughterTleTruncScan(
      output(), beamDaughterTleTruncDauIdx, beamDaughterTleTruncK, beamDaughterTleTruncNhitsInt,
      beamDaughterTleTruncTrueEkin0GeV, beamDaughterTleTruncPtleGeV, beamDaughterTleTruncFracRes);
  pionMomentumTree::FillPionMomentumVariables_BeamDaughterMcsTruncScan(
      output(), beamDaughterMcsTruncDauIdx, beamDaughterMcsTruncK, beamDaughterMcsTruncNsegments,
      beamDaughterMcsTruncTrueEkin0GeV, beamDaughterMcsTruncPmcsGeV, beamDaughterMcsTruncFracRes);

  if (_CreateEventDisplay && _EventDisplay) {
    _EventDisplay->FillTree(output(), GetEvent(), const_cast<ToyBoxPD*>(&box()));
  }
}

//********************************************************************
void pionMomentumAnalysis::FillToyVarsInMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
bool pionMomentumAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* /*part*/) {
//********************************************************************
  return false;
}

//********************************************************************
void pionMomentumAnalysis::FillTruthTree(const AnaTrueParticlePD& part) {
//********************************************************************
  pdBaseAnalysis::FillTruthTree(part);
}

//********************************************************************
void pionMomentumAnalysis::FillCategories() {
//********************************************************************
  if (box().MainTrack)
    anaUtils::FillCategories(&GetEvent(), box().MainTrack, "");
  pionMomentumAnaUtils::FillMainTrueStoppingCategory(box().MainTrack, _StoppingMaxTrueEndMomentumGeV);
}
