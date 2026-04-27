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
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {
bool BuildMCSLogLikelihoodFromSamples(const std::vector<double>& deltaTheta, const std::vector<double>& xOverX0,
                                      const std::vector<double>& pAxisGeV,
                                      const pdJointK0sPionMomentum::MCSLikelihoodConfig& cfg,
                                      std::vector<double>& logL) {
  logL.clear();
  if (deltaTheta.empty() || deltaTheta.size() != xOverX0.size() || pAxisGeV.empty()) return false;

  constexpr double kPionMassMeV = 139.57;
  constexpr double kHighlandConstantMeV = 13.6;
  const double theta0FloorRad =
      (std::isfinite(cfg.theta0FloorRad) && cfg.theta0FloorRad > 0.) ? cfg.theta0FloorRad : 1e-6;

  logL.reserve(pAxisGeV.size());
  for (double pGeV : pAxisGeV) {
    if (!std::isfinite(pGeV) || pGeV <= 0.0) return false;
    const double pMeV = pGeV * 1000.0;
    const double eMeV = std::sqrt(pMeV * pMeV + kPionMassMeV * kPionMassMeV);
    if (!std::isfinite(eMeV) || eMeV <= 0.0) return false;
    const double beta = pMeV / eMeV;
    if (!std::isfinite(beta) || beta <= 0.0) return false;

    double nll = 0.0;
    for (size_t i = 0; i < deltaTheta.size(); ++i) {
      const double xox0 = xOverX0[i];
      if (!std::isfinite(xox0) || xox0 <= 0.0) return false;
      double corr = 1.0 + 0.038 * std::log(std::max(xox0, 1e-12));
      if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;
      double theta0 = (kHighlandConstantMeV / (beta * pMeV)) * std::sqrt(xox0) * corr;
      if (!std::isfinite(theta0) || theta0 <= 0.0) theta0 = theta0FloorRad;
      theta0 = std::max(theta0, theta0FloorRad);
      const double pull = deltaTheta[i] / theta0;
      nll += 0.5 * pull * pull + std::log(theta0);
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
  _TLERequireTrueBeamPion = false;
  _EnsureMomentumSignalOnly = true;
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
  _TLERequireTrueBeamPion = ND::params().GetParameterI("pionMomentumAnalysis.TLERequireTrueBeamPion");
  _EnsureMomentumSignalOnly = ND::params().GetParameterI("pionMomentumAnalysis.EnsureMomentumSignalOnly");
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
  std::vector<double> mainsegmentlength;
  std::vector<double> mainthetascatterused;
  std::vector<double> mainrrscatterused;
  double mainmcsmomentum = -999.0;
  double maintlefitmomentum = -999.0;
  std::vector<double> pAxisTleGeV;
  std::vector<double> logLTle;
  std::vector<double> pAxisMcsGeV;
  std::vector<double> logLMcs;
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
  AnaParticlePD* beampart = NULL;
  if (beam && beam->BeamParticle) beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);

  if (box().MainTrack) {
    bool runTLE = true;
    if (_TLERequireTrueBeamPion) {
      runTLE = false;
      if (box().MainTrack->TrueObject && beampart && beampart->TrueObject) {
        const AnaTrueParticlePD* mainTrue = static_cast<const AnaTrueParticlePD*>(box().MainTrack->TrueObject);
        const AnaTrueParticlePD* beamTrue = static_cast<const AnaTrueParticlePD*>(beampart->TrueObject);
        if (mainTrue && beamTrue) {
          runTLE = (std::abs(mainTrue->PDG) == 211) && (mainTrue->ID == beamTrue->ID);
        }
      }
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

    pdJointK0sPionMomentum::MCSLikelihoodConfig mcsCfg;
    mcsCfg.radiationLengthCm = _MCSRadiationLengthCm;
    mcsCfg.targetSegmentLengthCm = _MCSTargetSegmentLengthCm;
    mcsCfg.minSegmentLengthCm = _MCSMinSegmentLengthCm;
    mcsCfg.theta0FloorRad = _MCStheta0FloorRad;
    mcsCfg.maxAbsDeltaThetaRad = _MCSMaxAbsDeltaThetaRad;

    std::vector<double> xoverx0;
    if (pdJointK0sPionMomentum::BuildPionMcsScatteringSamples(*box().MainTrack, mcsCfg, mainthetascatter, xoverx0, nullptr)) {
      mainsegmentlength.reserve(xoverx0.size());
      for (size_t i = 0; i < xoverx0.size(); ++i) {
        mainsegmentlength.push_back(xoverx0[i] * mcsCfg.radiationLengthCm);
      }

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
        mainrrscatterused.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        std::vector<double> xoverx0used;
        xoverx0used.reserve(static_cast<size_t>(lastKeepExcl - firstKeep));
        const double totalLen = std::accumulate(mainsegmentlength.begin(), mainsegmentlength.end(), 0.0);
        double accumLen = 0.0;
        for (int i = 0; i < firstKeep; ++i) accumLen += mainsegmentlength[static_cast<size_t>(i)];
        for (int i = firstKeep; i < lastKeepExcl; ++i) {
          const double seg = mainsegmentlength[static_cast<size_t>(i)];
          mainthetascatterused.push_back(mainthetascatter[static_cast<size_t>(i)]);
          const double segCenter = accumLen + 0.5 * seg;
          mainrrscatterused.push_back(std::max(0.0, totalLen - segCenter));
          accumLen += seg;
          xoverx0used.push_back(xoverx0[static_cast<size_t>(i)]);
        }

        pAxisMcsGeV.clear();
        for (double p = 0.05; p <= 2.50 + 1e-12; p += 0.01) pAxisMcsGeV.push_back(p);
        if (BuildMCSLogLikelihoodFromSamples(mainthetascatterused, xoverx0used, pAxisMcsGeV, mcsCfg, logLMcs) &&
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
      mainsegmentlength.clear();
      mainthetascatterused.clear();
      mainrrscatterused.clear();
      mainmcsmomentum = -999.0;
    }

    int stoppingCodeForDiag = -999;
    if (box().MainTrack->TrueObject) {
      const AnaTrueParticlePD* tpart = static_cast<const AnaTrueParticlePD*>(box().MainTrack->TrueObject);
      if (tpart) {
        if (std::abs(tpart->PDG) != 211) {
          stoppingCodeForDiag = 3;
        } else {
          const double pend = static_cast<double>(tpart->MomentumEnd);
          if (std::isfinite(pend) && pend >= 0.) {
            stoppingCodeForDiag = (pend <= _StoppingMaxTrueEndMomentumGeV) ? 1 : 2;
          }
        }
      }
    }
    pionMomentumAnaUtils::MomentumDiagConfig diagCfg;
    diagCfg.enableMomentumDiagnosticMultigraphs = _EnableMomentumDiagnosticMultigraphs;
    diagCfg.ensureMomentumSignalOnly = _EnsureMomentumSignalOnly;
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
    pionMomentumAnaUtils::MaybeAccumulateMainTrackMomentumDiagnostics(
        box().MainTrack, GetEvent(), stoppingCodeForDiag, diagCfg, maintlefitmomentum, mainmcsmomentum, pAxisTleGeV,
        logLTle, pAxisMcsGeV, logLMcs, mainthetascatter, mainsegmentlength, mainthetascatterused,
        mainrrscatterused);
  }

  pionMomentumTree::FillPionMomentumVariables_BeamParticleReco(output(), box().MainTrack, beampart,
                                                               mainthetascatter, mainsegmentlength, mainmcsmomentum,
                                                               maintlefitmomentum);

  if (_CreateEventDisplay && _EventDisplay) {
    _EventDisplay->FillTree(output(), GetEvent(), nullptr);
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
