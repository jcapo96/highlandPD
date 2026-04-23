#include "pionMomentumAnalysis.hxx"
#include "pionMomentumSelection.hxx"
#include "pionMomentumTree.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"
#include "pdJointK0sPionMomentum.hxx"
#include "pdUtilsTrack.hxx"

//********************************************************************
pionMomentumAnalysis::pionMomentumAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************
  _ApplySCECorrection = false;
  _ApplySCESystematic = false;
  _MCSRadiationLengthCm = 14.0;
  _MCSMinSegmentLengthCm = 0.5;
  _MCStheta0FloorRad = 1e-6;
  _MCSMaxAbsDeltaThetaRad = -1.0;
}

//********************************************************************
pionMomentumAnalysis::~pionMomentumAnalysis() {
//********************************************************************
}

//********************************************************************
bool pionMomentumAnalysis::Initialize() {
//********************************************************************
  if (!baseAnalysis::Initialize()) return false;

  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pionMomentumAnalysis.MinAccumLevelToSave"));
  _ApplySCECorrection = ND::params().GetParameterI("pionMomentumAnalysis.ApplySCECorrection");
  _ApplySCESystematic = ND::params().GetParameterI("pionMomentumAnalysis.ApplySCESystematic");
  _MCSRadiationLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSRadiationLengthCm");
  _MCSMinSegmentLengthCm = ND::params().GetParameterD("pionMomentumAnalysis.MCSMinSegmentLengthCm");
  _MCStheta0FloorRad = ND::params().GetParameterD("pionMomentumAnalysis.MCStheta0FloorRad");
  _MCSMaxAbsDeltaThetaRad = ND::params().GetParameterD("pionMomentumAnalysis.MCSMaxAbsDeltaThetaRad");

  return true;
}

//********************************************************************
void pionMomentumAnalysis::Finalize() {
//********************************************************************
}

//********************************************************************
void pionMomentumAnalysis::DefineInputConverters() {
//********************************************************************
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
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
  if (box().MainTrack) {
    // Ensure RR is available for hit-triplet ordering in MCS sampling.
    bool hasValidRR = false;
    for (size_t ihit = 0; ihit < box().MainTrack->Hits[2].size(); ++ihit) {
      const double rr = static_cast<double>(box().MainTrack->Hits[2][ihit].ResidualRange);
      if (std::isfinite(rr) && rr > 0.0) {
        hasValidRR = true;
        break;
      }
    }
    if (!hasValidRR && !box().MainTrack->Hits[2].empty()) pdAnaUtils::ComputeResidualRange(box().MainTrack);

    pdJointK0sPionMomentum::MCSLikelihoodConfig mcsCfg;
    mcsCfg.radiationLengthCm = _MCSRadiationLengthCm;
    mcsCfg.minSegmentLengthCm = _MCSMinSegmentLengthCm;
    mcsCfg.theta0FloorRad = _MCStheta0FloorRad;
    mcsCfg.maxAbsDeltaThetaRad = _MCSMaxAbsDeltaThetaRad;

    std::vector<double> xoverx0;
    if (pdJointK0sPionMomentum::BuildPionMcsScatteringSamples(*box().MainTrack, mcsCfg, mainthetascatter, xoverx0, nullptr)) {
      mainsegmentlength.reserve(xoverx0.size());
      for (size_t i = 0; i < xoverx0.size(); ++i) {
        mainsegmentlength.push_back(xoverx0[i] * mcsCfg.radiationLengthCm);
      }
    } else {
      mainthetascatter.clear();
      mainsegmentlength.clear();
    }
  }

  AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
  AnaParticlePD* beampart = NULL;
  if (beam && beam->BeamParticle) beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  pionMomentumTree::FillPionMomentumVariables_BeamParticleReco(output(), box().MainTrack, beampart,
                                                               mainthetascatter, mainsegmentlength);
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
}
