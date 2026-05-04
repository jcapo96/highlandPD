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
#include "pdJointK0sPionMomentum.hxx"
#include "pdUtilsDEdx.hxx"
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>

namespace {
bool BuildMCSLogLikelihoodFromSamples(const std::vector<double>& thetaXZ, const std::vector<double>& thetaYZ,
                                      const std::vector<double>& xOverX0, const std::vector<double>& pAxisGeV,
                                      const pdJointK0sPionMomentum::MCSLikelihoodConfig& cfg,
                                      std::vector<double>& logL) {
  logL.clear();
  if (thetaXZ.empty() || thetaYZ.size() != thetaXZ.size() || thetaXZ.size() != xOverX0.size() || pAxisGeV.empty())
    return false;

  constexpr double kPionMassMeV = 139.57;
  constexpr double kHighlandConstantMeV = 13.6;
  const double theta0FloorRad =
      (std::isfinite(cfg.theta0FloorRad) && cfg.theta0FloorRad > 0.) ? cfg.theta0FloorRad : 1e-6;
  const double sigmaDetFloorRad =
      (std::isfinite(cfg.detectorSigmaFloorRad) && cfg.detectorSigmaFloorRad > 0.) ? cfg.detectorSigmaFloorRad : 1e-6;
  const double sigmaDetA = (std::isfinite(cfg.detectorSigmaA) && cfg.detectorSigmaA >= 0.) ? cfg.detectorSigmaA : 0.0;
  const double sigmaDetC = std::isfinite(cfg.detectorSigmaC) ? cfg.detectorSigmaC : 0.0;
  const double x0Cm = (std::isfinite(cfg.radiationLengthCm) && cfg.radiationLengthCm > 1e-9) ? cfg.radiationLengthCm : 14.0;

  logL.reserve(pAxisGeV.size());
  for (double pGeV : pAxisGeV) {
    if (!std::isfinite(pGeV) || pGeV <= 0.0) return false;
    const double pMeV = pGeV * 1000.0;
    const double eMeV = std::sqrt(pMeV * pMeV + kPionMassMeV * kPionMassMeV);
    if (!std::isfinite(eMeV) || eMeV <= 0.0) return false;
    const double beta = pMeV / eMeV;
    if (!std::isfinite(beta) || beta <= 0.0) return false;

    double nll = 0.0;
    for (size_t i = 0; i < thetaXZ.size(); ++i) {
      const double xox0 = xOverX0[i];
      if (!std::isfinite(xox0) || xox0 <= 0.0) return false;
      double corr = 1.0 + 0.038 * std::log(std::max(xox0, 1e-12));
      if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;
      double theta0 = (kHighlandConstantMeV / (beta * pMeV)) * std::sqrt(xox0) * corr;
      if (!std::isfinite(theta0) || theta0 <= 0.0) theta0 = theta0FloorRad;
      theta0 = std::max(theta0, theta0FloorRad);
      double sigma = theta0;
      if (cfg.useDetectorSigma) {
        const double segCm = xox0 * x0Cm;
        double sigmaDet2 = 0.0;
        if (std::isfinite(segCm) && segCm > 0.0) sigmaDet2 = sigmaDetA / (segCm * segCm) + sigmaDetC;
        if (!std::isfinite(sigmaDet2) || sigmaDet2 <= 0.0) sigmaDet2 = sigmaDetFloorRad * sigmaDetFloorRad;
        const double sigmaDet = std::max(std::sqrt(sigmaDet2), sigmaDetFloorRad);
        sigma = std::sqrt(theta0 * theta0 + sigmaDet * sigmaDet);
      }
      sigma = std::max(sigma, theta0FloorRad);
      const double dtxz = thetaXZ[i];
      const double dtyz = thetaYZ[i];
      if (!std::isfinite(dtxz) || !std::isfinite(dtyz)) return false;
      const double pullX = dtxz / sigma;
      const double pullY = dtyz / sigma;
      nll += 0.5 * (pullX * pullX + pullY * pullY) + 2.0 * std::log(sigma);
    }
    if (!std::isfinite(nll)) return false;
    logL.push_back(-nll);
  }
  return logL.size() == pAxisGeV.size();
}
}  // namespace

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
  _BraggChi2PiSigmaMeVcm = 0.3;
  _BraggChi2PiMinHits = 3;
  _BraggChi2PiSkipHitsFirst = 0;
  _BraggChi2PiSkipHitsLast = 0;
  _BraggChi2PiDedxMinMeVcm = 0.;
  _BraggChi2PiDedxMaxMeVcm = 0.;
  _MomentumRequireTrueBeamPion = false;
  _RunTLEFit = true;
  _RunMCS = true;
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
  _MCSRadiationLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSRadiationLengthCm");
  _MCSTargetSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSTargetSegmentLengthCm");
  _MCSMinSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSMinSegmentLengthCm");
  _MCStheta0FloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCStheta0FloorRad");
  _MCSMaxAbsDeltaThetaRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSMaxAbsDeltaThetaRad");
  _MCSDropFirstNSegments = ND::params().GetParameterI("pionMomentumAnalysis.MCSDropFirstNSegments");
  _MCSDropLastNSegments = ND::params().GetParameterI("pionMomentumAnalysis.MCSDropLastNSegments");
  _MCSUseDetectorSigma = ND::params().GetParameterI("pionMomentumAnalysis.MCSUseDetectorSigma") != 0;
  _MCSDetectorSigmaFloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaFloorRad");
  _MCSDetectorSigmaA = ND::params().HasParameter("pionMomentumAnalysis.MCSDetectorSigmaA")
                           ? ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaA")
                           : 0.0;
  _MCSDetectorSigmaC = ND::params().HasParameter("pionMomentumAnalysis.MCSDetectorSigmaC")
                           ? ND::params().GetParameterD("pionMomentumAnalysis.MCSDetectorSigmaC")
                           : 0.0;
  _MCSMomentumScanMaxGeV = ND::params().HasParameter("pionMomentumAnalysis.MCSMomentumScanMaxGeV")
                               ? ND::params().GetParameterD("pionMomentumAnalysis.MCSMomentumScanMaxGeV")
                               : 2.5;
  _TLEMinInteriorHits = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxMinInteriorHits");
  _TLESkipHitsFirst = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxSkipHitsFirst");
  _TLESkipHitsLast = ND::params().GetParameterI("pionMomentumAnalysis.FreeRangeDedxSkipHitsLast");
  _TLEDedxMinMeVcm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxMinMeVcm");
  _TLEDedxMaxMeVcm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxMaxMeVcm");
  _TLEScanLmaxCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanLmaxCm");
  _TLEScanStepCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanStepCm");
  _TLEScanStepFineCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeScanStepFineCm");
  _TLELowPMomentumRefineGeV = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeLowPMomentumRefineGeV");
  _TLEDedxPdfPathCm = ND::params().GetParameterD("pionMomentumAnalysis.FreeRangeDedxPdfPathCm");
  _BraggChi2PiMaxResidualRangeCm = ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiMaxResidualRangeCm");
  _BraggChi2PiSigmaMeVcm = ND::params().GetParameterD("pionMomentumAnalysis.BraggChi2PiSigmaMeVcm");
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

    std::vector<double> rrMain;
    std::vector<double> dedxMain;
    if (runTLE) {
      rrMain.reserve(box().MainTrack->Hits[2].size());
      dedxMain.reserve(box().MainTrack->Hits[2].size());
      for (size_t ihit = 0; ihit < box().MainTrack->Hits[2].size(); ++ihit) {
        const double rr = static_cast<double>(box().MainTrack->Hits[2][ihit].ResidualRange);
        const double dedx = static_cast<double>(box().MainTrack->Hits[2][ihit].dEdx);
        if (!std::isfinite(rr) || !std::isfinite(dedx)) continue;
        rrMain.push_back(rr);
        dedxMain.push_back(dedx);
      }

      if (pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurveFromVectors(
              rrMain, dedxMain, static_cast<double>(box().MainTrack->Length), _TLEScanLmaxCm, _TLEScanStepCm,
              _TLEMinInteriorHits, _TLESkipHitsFirst, _TLESkipHitsLast, _TLEDedxMinMeVcm, _TLEDedxMaxMeVcm,
              _TLEDedxPdfPathCm, pAxisTleGeV, logLTle, _TLEScanStepFineCm, _TLELowPMomentumRefineGeV) &&
          logLTle.size() == pAxisTleGeV.size()) {
        double bestLogL = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < pAxisTleGeV.size(); ++i) {
          if (!std::isfinite(logLTle[i])) continue;
          if (logLTle[i] > bestLogL) {
            bestLogL = logLTle[i];
            maintlefitmomentum = pAxisTleGeV[i];
          }
        }
      }
    }

    if (runMCS) {
    pdJointK0sPionMomentum::MCSLikelihoodConfig mcsCfg;
    mcsCfg.radiationLengthCm = _MCSRadiationLengthCm;
    mcsCfg.targetSegmentLengthCm = _MCSTargetSegmentLengthCm;
    mcsCfg.minSegmentLengthCm = _MCSMinSegmentLengthCm;
    mcsCfg.theta0FloorRad = _MCStheta0FloorRad;
    mcsCfg.maxAbsDeltaThetaRad = _MCSMaxAbsDeltaThetaRad;
    mcsCfg.useDetectorSigma = _MCSUseDetectorSigma;
    mcsCfg.detectorSigmaFloorRad = _MCSDetectorSigmaFloorRad;
    mcsCfg.detectorSigmaA = _MCSDetectorSigmaA;
    mcsCfg.detectorSigmaC = _MCSDetectorSigmaC;

    std::vector<double> xoverx0;
    if (pdJointK0sPionMomentum::BuildPionMcsScatteringSamples(*box().MainTrack, mcsCfg, mainthetascatter, xoverx0,
                                                              nullptr, &mainthetascatterxz, &mainthetascatteryz,
                                                              nullptr, nullptr, nullptr, &mainsegmentlengthsingle,
                                                              &mainsegmentfitchi2ndfsingle, &mainsegmentrrtoendsingle)) {
      mainsegmentlengthpair.clear();
      mainsegmentlengthpair.reserve(xoverx0.size());
      for (size_t i = 0; i < xoverx0.size(); ++i) mainsegmentlengthpair.push_back(xoverx0[i] * mcsCfg.radiationLengthCm);

      const int nMcs = static_cast<int>(mainthetascatter.size());
      int dropFirst = 0;
      int dropLast = 0;
      if (nMcs > 0) {
        dropFirst = std::max(0, _MCSDropFirstNSegments);
        dropLast = std::max(0, _MCSDropLastNSegments);
        dropFirst = std::min(dropFirst, nMcs - 1);
        dropLast = std::min(dropLast, nMcs - 1);
        if (dropFirst + dropLast >= nMcs) {
          dropFirst = 0;
          dropLast = 0;
        }
      }
      const int firstKeep = dropFirst;
      const int lastKeepExcl = nMcs - dropLast;
      if (nMcs > 0 && lastKeepExcl > firstKeep) {
        mainthetascatterused.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        mainthetascatterxzusd.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        mainthetascatteryzusd.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        mainrrscatterused.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        std::vector<double> xoverx0used;
        xoverx0used.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        const double totalLen = std::accumulate(mainsegmentlengthpair.begin(), mainsegmentlengthpair.end(), 0.0);
        double accumLen = 0.0;
        for (int i = 0; i < firstKeep; ++i) accumLen += mainsegmentlengthpair[static_cast<size_t>(i)];
        for (int i = firstKeep; i < lastKeepExcl; ++i) {
          const double seg = mainsegmentlengthpair[static_cast<size_t>(i)];
          mainthetascatterused.push_back(mainthetascatter[static_cast<size_t>(i)]);
          mainthetascatterxzusd.push_back(mainthetascatterxz[static_cast<size_t>(i)]);
          mainthetascatteryzusd.push_back(mainthetascatteryz[static_cast<size_t>(i)]);
          const double segCenter = accumLen + 0.5 * seg;
          mainrrscatterused.push_back(std::max(0.0, totalLen - segCenter));
          accumLen += seg;
          xoverx0used.push_back(xoverx0[static_cast<size_t>(i)]);
        }

        pAxisMcsGeV.clear();
        const double pMcsMax = (std::isfinite(_MCSMomentumScanMaxGeV) && _MCSMomentumScanMaxGeV > 0.05)
                                   ? _MCSMomentumScanMaxGeV
                                   : 2.5;
        for (double p = 0.05; p <= pMcsMax + 1e-12; p += 0.01) pAxisMcsGeV.push_back(p);
        if (BuildMCSLogLikelihoodFromSamples(mainthetascatterxzusd, mainthetascatteryzusd, xoverx0used, pAxisMcsGeV,
                                             mcsCfg, logLMcs) &&
            logLMcs.size() == pAxisMcsGeV.size()) {
          double bestLogL = -std::numeric_limits<double>::infinity();
          for (size_t i = 0; i < pAxisMcsGeV.size(); ++i) {
            if (!std::isfinite(logLMcs[i])) continue;
            if (logLMcs[i] > bestLogL) {
              bestLogL = logLMcs[i];
              mainmcsmomentum = pAxisMcsGeV[i];
            }
          }
        }
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
    std::vector<double> maintruexoverx0;
    if (pdJointK0sPionMomentum::BuildPionMcsScatteringSamplesFromOrderedPositions(
            trueOrdered, mcsCfg, maintruethetascatter, maintruexoverx0, nullptr, &maintruethetascatterxz,
            &maintruethetascatteryz, nullptr, nullptr, nullptr, &maintruesegmentlengthsingle,
            &maintruesegmentfitchi2ndfsingle, &maintruesegmentrrtoendsingle)) {
      maintruesegmentlengthpair.clear();
      maintruesegmentlengthpair.reserve(maintruexoverx0.size());
      for (size_t i = 0; i < maintruexoverx0.size(); ++i)
        maintruesegmentlengthpair.push_back(maintruexoverx0[i] * mcsCfg.radiationLengthCm);

      const int nMcsT = static_cast<int>(maintruethetascatter.size());
      int dropFirstT = 0;
      int dropLastT = 0;
      if (nMcsT > 0) {
        dropFirstT = std::max(0, _MCSDropFirstNSegments);
        dropLastT = std::max(0, _MCSDropLastNSegments);
        dropFirstT = std::min(dropFirstT, nMcsT - 1);
        dropLastT = std::min(dropLastT, nMcsT - 1);
        if (dropFirstT + dropLastT >= nMcsT) {
          dropFirstT = 0;
          dropLastT = 0;
        }
      }
      const int firstKeepT = dropFirstT;
      const int lastKeepExclT = nMcsT - dropLastT;
      if (nMcsT > 0 && lastKeepExclT > firstKeepT) {
        std::vector<double> maintruethetascatterxzusd;
        std::vector<double> maintruethetascatteryzusd;
        std::vector<double> xoverx0usedT;
        maintruethetascatterxzusd.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
        maintruethetascatteryzusd.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
        xoverx0usedT.reserve(static_cast<size_t>(lastKeepExclT - firstKeepT));
        for (int i = firstKeepT; i < lastKeepExclT; ++i) {
          maintruethetascatterxzusd.push_back(maintruethetascatterxz[static_cast<size_t>(i)]);
          maintruethetascatteryzusd.push_back(maintruethetascatteryz[static_cast<size_t>(i)]);
          xoverx0usedT.push_back(maintruexoverx0[static_cast<size_t>(i)]);
        }
        std::vector<double> pAxisMcsTrueGeV;
        std::vector<double> logLMcsTrue;
        const double pMcsMaxT = (std::isfinite(_MCSMomentumScanMaxGeV) && _MCSMomentumScanMaxGeV > 0.05)
                                    ? _MCSMomentumScanMaxGeV
                                    : 2.5;
        for (double p = 0.05; p <= pMcsMaxT + 1e-12; p += 0.01) pAxisMcsTrueGeV.push_back(p);
        if (BuildMCSLogLikelihoodFromSamples(maintruethetascatterxzusd, maintruethetascatteryzusd, xoverx0usedT,
                                             pAxisMcsTrueGeV, mcsCfg, logLMcsTrue) &&
            logLMcsTrue.size() == pAxisMcsTrueGeV.size()) {
          double bestLogLT = -std::numeric_limits<double>::infinity();
          for (size_t i = 0; i < pAxisMcsTrueGeV.size(); ++i) {
            if (!std::isfinite(logLMcsTrue[i])) continue;
            if (logLMcsTrue[i] > bestLogLT) {
              bestLogLT = logLMcsTrue[i];
              maintruemcsmomentum = pAxisMcsTrueGeV[i];
            }
          }
        }
      }
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
    if (_EnableMomentumDiagnosticMultigraphs && (runTLE || runMCS)) {
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

  // Bragg χ² (arXiv:2409.18288 Eq. 6.1): runs whenever MainTrack exists — not gated on RunTLEFit or RunMCS.
  double mainbraggchi2pibb = -999.;
  Int_t mainbraggdedxnhits = -999;
  if (box().MainTrack) {
    double chi2 = -999.;
    int nh = 0;
    if (pdAnaUtils::ComputePionBraggWindowChi2PiEq61(box().MainTrack, _BraggChi2PiMaxResidualRangeCm,
                                                     _BraggChi2PiSigmaMeVcm, _BraggChi2PiMinHits,
                                                     _BraggChi2PiSkipHitsFirst, _BraggChi2PiSkipHitsLast,
                                                     _BraggChi2PiDedxMinMeVcm, _BraggChi2PiDedxMaxMeVcm, chi2, nh)) {
      mainbraggchi2pibb = chi2;
      mainbraggdedxnhits = static_cast<Int_t>(nh);
    }
  }

  pionMomentumTree::FillPionMomentumVariables_BeamParticleReco(
      output(), box().MainTrack, beampart, mainthetascatter, mainthetascatterxz, mainthetascatteryz,
      mainsegmentlengthsingle, mainsegmentfitchi2ndfsingle, mainsegmentrrtoendsingle, mainmcsmomentum,
      maintlefitmomentum, mainbraggchi2pibb, mainbraggdedxnhits, mainntruetrajpointsTree, maintruethetascatter,
      maintruethetascatterxz, maintruethetascatteryz, maintruesegmentlengthsingle, maintruesegmentfitchi2ndfsingle,
      maintruesegmentrrtoendsingle, maintruemcsmomentum);

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
