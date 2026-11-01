#include "neutralKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "pdMomReconstructionFromParams.hxx"
#include "neutralKaonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"
#include "neutralKaonTruthTree.hxx"
#include "ToyBoxNeutralKaon.hxx"

#include "pdAnalysisUtils.hxx"
#include "pdUtilsRangePID.hxx"
#include "pdAnnihilationUtils.hxx"
#include "standardPDTree.hxx"
#include <sstream>
#include <map>
#include <set>

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

#include "CategoryManager.hxx"

#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"

#include "baseToyMaker.hxx"

#include "AnalysisUtils.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include "TVector3.h"

namespace {
  bool HasValidRecoStart(const AnaParticlePD* recoPart){
    if(!recoPart) return false;
    return recoPart->PositionStart[0] > -900.f &&
           recoPart->PositionStart[1] > -900.f &&
           recoPart->PositionStart[2] > -900.f;
  }

  Float_t RecoDaughtersStartDistance(const AnaParticlePD* r1, const AnaParticlePD* r2){
    if(!HasValidRecoStart(r1) || !HasValidRecoStart(r2)) return static_cast<Float_t>(-999.0);
    const double dx = static_cast<double>(r1->PositionStart[0] - r2->PositionStart[0]);
    const double dy = static_cast<double>(r1->PositionStart[1] - r2->PositionStart[1]);
    const double dz = static_cast<double>(r1->PositionStart[2] - r2->PositionStart[2]);
    return static_cast<Float_t>(std::sqrt(dx * dx + dy * dy + dz * dz));
  }

  bool HasValidTrueStart(const AnaTrueParticlePD* truePart){
    if(!truePart) return false;
    return truePart->Position[0] > -900.f &&
           truePart->Position[1] > -900.f &&
           truePart->Position[2] > -900.f;
  }

  bool HasValidTrueEnd(const AnaTrueParticlePD* truePart){
    if(!truePart) return false;
    return truePart->PositionEnd[0] > -900.f &&
           truePart->PositionEnd[1] > -900.f &&
           truePart->PositionEnd[2] > -900.f;
  }

  Float_t TrueScalarMomentum(const AnaTrueParticlePD* truePart, bool atEnd){
    if(!truePart) return static_cast<Float_t>(-999.0);
    const Float_t m = atEnd ? truePart->MomentumEnd : truePart->Momentum;
    if(!std::isfinite(m)) return static_cast<Float_t>(-999.0);
    return m;
  }

  // Rest mass (MeV) for common MC secondaries; unknown PDG -> 0 (E ~ |p| for highly relativistic).
  double TrueRestMassMeV(Int_t pdg){
    const int ap = std::abs(static_cast<int>(pdg));
    switch(ap){
      case 211:
        return 139.57;
      case 111:
        return 134.98;
      case 2212:
        return 938.27;
      case 2112:
        return 939.57;
      case 11:
        return 0.511;
      case 13:
        return 105.66;
      case 321:
        return 493.68;
      case 310:
      case 130:
        return 497.61;
      case 22:
        return 0.0;
      case 12:
      case 14:
      case 16:
        return 0.0;
      default:
        return 0.0;
    }
  }

  Float_t TrueTotalEnergyGeV(const AnaTrueParticlePD* truePart, bool atEnd){
    if(!truePart) return static_cast<Float_t>(-999.0);
    const double p = static_cast<double>(atEnd ? truePart->MomentumEnd : truePart->Momentum);
    if(!std::isfinite(p) || p < 0.0) return static_cast<Float_t>(-999.0);
    const double mGeV = TrueRestMassMeV(truePart->PDG) / 1000.0;
    return static_cast<Float_t>(std::sqrt(p * p + mGeV * mGeV));
  }

  // Same as neutralKaonTree: |PositionEnd - Position| when both vertices are valid (not AnaTrueParticle::Length).
  bool IsValidTruePos3(Float_t x, Float_t y, Float_t z){
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z) && x > -900.f && y > -900.f && z > -900.f;
  }

  Float_t TrueLengthVertexSeparationCm(const AnaTrueParticlePD* truePart){
    if(!truePart) return static_cast<Float_t>(-999.0);
    const TVector3 trueStart(truePart->Position[0], truePart->Position[1], truePart->Position[2]);
    const TVector3 trueEnd(truePart->PositionEnd[0], truePart->PositionEnd[1], truePart->PositionEnd[2]);
    if(!IsValidTruePos3(trueStart.X(), trueStart.Y(), trueStart.Z()) ||
       !IsValidTruePos3(trueEnd.X(), trueEnd.Y(), trueEnd.Z()))
      return static_cast<Float_t>(-999.0);
    return static_cast<Float_t>((trueEnd - trueStart).Mag());
  }

  Float_t RecoLengthOrSentinel(const AnaParticlePD* recoPart){
    if(!recoPart) return static_cast<Float_t>(-999.0);
    const Float_t L = recoPart->Length;
    if(!std::isfinite(L) || L < 0.f) return static_cast<Float_t>(-999.0);
    return L;
  }

  Int_t RecoNHitsOrSentinel(const AnaParticlePD* recoPart){
    if(!recoPart) return -999;
    if(recoPart->NHits < 0) return -999;
    return recoPart->NHits;
  }

  Int_t TruePdgOfReco(const AnaParticlePD* recoPart){
    if(!recoPart || !recoPart->TrueObject) return -999;
    auto* t = static_cast<AnaTrueParticlePD*>(recoPart->TrueObject);
    return t ? t->PDG : -999;
  }

  Float_t RecoVsTrueStartDistance(const AnaParticlePD* recoPart, const AnaTrueParticlePD* truePart){
    if(!HasValidRecoStart(recoPart) || !HasValidTrueStart(truePart)) return static_cast<Float_t>(-999.0);
    const double dx = static_cast<double>(recoPart->PositionStart[0] - truePart->Position[0]);
    const double dy = static_cast<double>(recoPart->PositionStart[1] - truePart->Position[1]);
    const double dz = static_cast<double>(recoPart->PositionStart[2] - truePart->Position[2]);
    return static_cast<Float_t>(std::sqrt(dx * dx + dy * dy + dz * dz));
  }

  Float_t RecoVsTrueEndDistance(const AnaParticlePD* recoPart, const AnaTrueParticlePD* truePart){
    if(!recoPart || !truePart) return static_cast<Float_t>(-999.0);
    const bool validRecoEnd = (recoPart->PositionEnd[0] > -900.f &&
                               recoPart->PositionEnd[1] > -900.f &&
                               recoPart->PositionEnd[2] > -900.f);
    if(!validRecoEnd || !HasValidTrueEnd(truePart)) return static_cast<Float_t>(-999.0);
    const double dx = static_cast<double>(recoPart->PositionEnd[0] - truePart->PositionEnd[0]);
    const double dy = static_cast<double>(recoPart->PositionEnd[1] - truePart->PositionEnd[1]);
    const double dz = static_cast<double>(recoPart->PositionEnd[2] - truePart->PositionEnd[2]);
    return static_cast<Float_t>(std::sqrt(dx * dx + dy * dy + dz * dz));
  }

  bool IsValidDistanceValue(Float_t value){
    return std::isfinite(value) && value > -900.f;
  }

  Int_t VertexGlobalMinIsDauPairTag(Float_t dauPair, Float_t d1Other, Float_t d2Other){
    if(!IsValidDistanceValue(dauPair)) return 0;
    Float_t globalMin = dauPair;
    if(IsValidDistanceValue(d1Other)) globalMin = std::min(globalMin, d1Other);
    if(IsValidDistanceValue(d2Other)) globalMin = std::min(globalMin, d2Other);
    return (std::abs(dauPair - globalMin) <= 1e-5f) ? 1 : 0;
  }

  Float_t RecoPairFitMinimumDistance(AnaParticlePD* p1, AnaParticlePD* p2,
                                     double trackFitLength, double trackFitDistanceFromStart){
    if(!HasValidRecoStart(p1) || !HasValidRecoStart(p2)) return static_cast<Float_t>(-999.0);

    std::vector<double> fit1;
    std::vector<double> fit2;
    pdAnaUtils::ExtrapolateTrack(p1, fit1, trackFitLength, true, trackFitDistanceFromStart);
    pdAnaUtils::ExtrapolateTrack(p2, fit2, trackFitLength, true, trackFitDistanceFromStart);

    const bool fit1Valid =
      (fit1.size() >= 6 &&
       std::isfinite(fit1[0]) && std::isfinite(fit1[1]) && std::isfinite(fit1[2]) &&
       std::isfinite(fit1[3]) && std::isfinite(fit1[4]) && std::isfinite(fit1[5]) &&
       fit1[0] > -900.0);
    const bool fit2Valid =
      (fit2.size() >= 6 &&
       std::isfinite(fit2[0]) && std::isfinite(fit2[1]) && std::isfinite(fit2[2]) &&
       std::isfinite(fit2[3]) && std::isfinite(fit2[4]) && std::isfinite(fit2[5]) &&
       fit2[0] > -900.0);
    if(!fit1Valid || !fit2Valid) return static_cast<Float_t>(-999.0);

    TVector3 closest1;
    TVector3 closest2;
    const double minDistance = pdAnaUtils::FindClosestPointsBetweenLines(fit1, fit2, closest1, closest2);
    if(!std::isfinite(minDistance)) return static_cast<Float_t>(-999.0);
    return static_cast<Float_t>(minDistance);
  }

  Float_t RecoProtonChi2Ndf(const AnaParticlePD* recoPart){
    return pdAnaUtils::Chi2PIDChi2PerHit(recoPart, 2212);
  }

  Float_t RecoKaonChi2Ndf(const AnaParticlePD* recoPart){
    return pdAnaUtils::Chi2PIDChi2PerHit(recoPart, 321);
  }

  Float_t RecoPionChi2Ndf(const AnaParticlePD* recoPart){
    return pdAnaUtils::Chi2PIDChi2PerHit(recoPart, 211);
  }

  Int_t RecoDaughterParentUidMinusTrueParentRecoUid(const AnaBunchB& bunch,
                                                    const AnaParticlePD* daughterReco,
                                                    const AnaParticlePD* trueParentReco){
    if(!daughterReco || !trueParentReco) return -999;
    if(daughterReco->ParentID < 0) return -999;
    AnaParticleB* parentReco = anaUtils::GetParticleByID(bunch, daughterReco->ParentID);
    if(!parentReco) return -999;
    return parentReco->UniqueID - trueParentReco->UniqueID;
  }

  bool HasValidRecoPositions(const AnaParticlePD* recoPart){
    if(!recoPart) return false;

    const bool validStart = (recoPart->PositionStart[0] > -900.f &&
                             recoPart->PositionStart[1] > -900.f &&
                             recoPart->PositionStart[2] > -900.f);
    const bool validEnd = (recoPart->PositionEnd[0] > -900.f &&
                           recoPart->PositionEnd[1] > -900.f &&
                           recoPart->PositionEnd[2] > -900.f);
    return validStart && validEnd;
  }

  Int_t RecoZStartGreaterThanEndFlag(const AnaParticlePD* recoPart){
    if(!recoPart) return 0;
    const Float_t zs = recoPart->PositionStart[2];
    const Float_t ze = recoPart->PositionEnd[2];
    if(zs <= -900.f || ze <= -900.f) return 0;
    return (zs > ze) ? 1 : 0;
  }

  bool HasValidRecoHits(const AnaParticlePD* recoPart){
    if(!recoPart) return false;

    const int nHitsTotal = recoPart->NHitsPerPlane[0] + recoPart->NHitsPerPlane[1] + recoPart->NHitsPerPlane[2];
    if(nHitsTotal > 0) return true;

    return (!recoPart->Hits[0].empty() || !recoPart->Hits[1].empty() || !recoPart->Hits[2].empty());
  }

  bool HasRecoObjectForTruthFlags(const AnaParticlePD* recoPart){
    return HasValidRecoPositions(recoPart) && HasValidRecoHits(recoPart);
  }
}

//********************************************************************
neutralKaonAnalysis::neutralKaonAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************

  _EventDisplayDataTree = nullptr;

  // Add the package version
  //  ND::versioning().AddPackage("StoppingProtonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
}

//********************************************************************
neutralKaonAnalysis::~neutralKaonAnalysis() {
//********************************************************************
  if (_EventDisplayDataTree) {
    delete _EventDisplayDataTree;
    _EventDisplayDataTree = nullptr;
  }
}

//********************************************************************
bool neutralKaonAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file.
     This is because in order to the overwride parameters file
     option (-p param.dat) to work, parameters cannot be accessed in the constructors.
  */

  // Initialize the pdBaseAnalysis
  if(!pdBaseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("neutralKaonAnalysis.MinAccumLevelToSave"));

  // SCE correction parameter
  _ApplySCECorrection = ND::params().GetParameterI("neutralKaonAnalysis.ApplySCECorrection");
  _ApplySCESystematic = ND::params().GetParameterI("neutralKaonAnalysis.ApplySCESystematic");

  // dE/dx likelihood model of the TLE fit (off by default: legacy PDF)
  pdMomReconstruction::SetDedxPdfModel_FromParams("neutralKaonAnalysis");

  // Initialize event display data tree
  // This stores data for later visualization via DrawingTools in ROOT macros
  _EventDisplayDataTree = new EventDisplayDataTree();

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("bestcandidate");

  // Add custom categories for neutral particle analysis
  neutralKaonAnaUtils::AddCustomCategories();

  return true;
}

//********************************************************************
void neutralKaonAnalysis::Finalize(){
//********************************************************************
  // Persist external diagnostic profiles filled during event loop.
  neutralKaonTree::WriteHitDistanceProfiles(output());
}

//********************************************************************
void neutralKaonAnalysis::DefineInputConverters(){
//********************************************************************

  /* In this method we add the to the InputManager (accessed by input() ) the InputConverters created
     in separet files (see, for example, pdIO/PDSPAnalyzerTreeConverter.hxx)
     which define the allowed input file formats.
  */

  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void neutralKaonAnalysis::DefineSelections(){
//********************************************************************

  sel().AddSelection("neutralKaonSelection", "Neutral Kaon Selection", new neutralKaonSelection(false)); // true/false for forcing break
}

//********************************************************************
void neutralKaonAnalysis::DefineCorrections(){
//********************************************************************

  // Some corrections are defined in pdBaseAnalysis
  pdBaseAnalysis::DefineCorrections();
  if(_ApplySCECorrection)
    corr().AddCorrection(0, "sce geometric correction", new ParticlePositionSCECorrection());
}

//********************************************************************
void neutralKaonAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in pdBaseAnalysis (highland/src/highland2/pdBaseAnalysis)
  pdBaseAnalysis::DefineSystematics();

  if(_ApplySCESystematic)
    evar().AddEventVariation(kSCEGeometric, "SCE variation", new SCEGeometricVariation());
}

//********************************************************************
void neutralKaonAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in pdBaseAnalysis
  pdBaseAnalysis::DefineConfigurations();
}

//********************************************************************
void neutralKaonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  if (addBase) pdBaseAnalysis::DefineMicroTrees(addBase);

  // // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_EventInfo(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  // standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  // standardPDTree::AddStandardVariables_BeamParticleReco(output());
  // standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
  // standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),50);
  // standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),50);
  // standardPDTree::AddStandardVariables_BeamTruthDaughters(output(),50);

  // AddVarI(output(), seltrk_dau_trueparentpdg, "Parent PDG of reco daughter");

  // AddVarI(output(), nAllParticles, "Number of all particles with valid start positions");

  // AddVarI would name the branches after the qualified enum ("neutralKaonTree::..."), hence AddVar
  output().AddVar(neutralKaonTree::ntruelambdapairs, "ntruelambdapairs", "I",
                  "True p + pi- with the same parent, both created by a decay (Lambda -> p pi-), among stored true "
                  "particles", -9999);
  output().AddVar(neutralKaonTree::ntruelambdapairsreco, "ntruelambdapairsreco", "I",
                  "Those pairs with both daughters matched to a reco object", -9999);
  output().AddVar(neutralKaonTree::lampairrecoflag, "lampairrecoflag", "I",
                  "Lambda pair: -1 none, 0 both daughters reco, 1 proton not reco, 2 pion not reco, 3 neither", -9999);
  output().AddVar(neutralKaonTree::lampairstage, "lampairstage", "I",
                  "Lambda pair: first failed vertex requirement (0 none, 1 hits, 2 parent, 3 PID, 4 distance, 5 overlap)",
                  -9999);
  output().AddVar(neutralKaonTree::lampairfailmask, "lampairfailmask", "I",
                  "Lambda pair: every failed requirement, bits 1 hits, 2 parent, 4 PID, 8 distance", -9999);
  output().AddVar(neutralKaonTree::lampairdist, "lampairdist", "F", "Lambda pair: smallest endpoint distance [cm]",
                  -9999);
  output().AddVar(neutralKaonTree::lampairnhitsp, "lampairnhitsp", "I", "Lambda pair: proton collection hits", -9999);
  output().AddVar(neutralKaonTree::lampairnhitspi, "lampairnhitspi", "I", "Lambda pair: pion collection hits", -9999);
  output().AddVar(neutralKaonTree::lampairpchi2p, "lampairpchi2p", "F", "Lambda pair: proton Chi2PID(2212)/npts",
                  -9999);
  output().AddVar(neutralKaonTree::lampairpichi2p, "lampairpichi2p", "F", "Lambda pair: pion Chi2PID(2212)/npts",
                  -9999);
  output().AddVar(neutralKaonTree::lampairpichi2pi, "lampairpichi2pi", "F", "Lambda pair: pion Chi2PID(211)/npts",
                  -9999);
  output().AddVar(neutralKaonTree::lampairtruemomp, "lampairtruemomp", "F", "Lambda pair: true proton momentum",
                  -9999);
  output().AddVar(neutralKaonTree::lampairtruemompi, "lampairtruemompi", "F", "Lambda pair: true pion momentum",
                  -9999);
  output().AddVar(neutralKaonTree::lampairtruelengthp, "lampairtruelengthp", "F",
                  "Lambda pair: true proton start-end distance [cm]", -9999);
  output().AddVar(neutralKaonTree::lampairtruelengthpi, "lampairtruelengthpi", "F",
                  "Lambda pair: true pion start-end distance [cm]", -9999);
  output().AddVar(neutralKaonTree::lampairptrueendproc, "lampairptrueendproc", "I",
                  "Lambda pair: true proton ProcessEnd (8 protonInelastic)", -9999);
  output().AddVar(neutralKaonTree::lampairpitrueendproc, "lampairpitrueendproc", "I",
                  "Lambda pair: true pion ProcessEnd (2 Decay, 10 piminusInelastic, 11 capture at rest)", -9999);
  output().AddVar(neutralKaonTree::lampairprecoreversed, "lampairprecoreversed", "I",
                  "Lambda pair: reco proton runs against the true one (1), along it (0)", -9999);
  output().AddVar(neutralKaonTree::lampairpirecoreversed, "lampairpirecoreversed", "I",
                  "Lambda pair: reco pion runs against the true one (1), along it (0)", -9999);
  output().AddVar(neutralKaonTree::lampairpvtxreversed, "lampairpvtxreversed", "I",
                  "Lambda pair: the proton endpoint nearest the pion is its reco end (1) or start (0)", -9999);
  output().AddVar(neutralKaonTree::lampairpivtxreversed, "lampairpivtxreversed", "I",
                  "Lambda pair: the pion endpoint nearest the proton is its reco end (1) or start (0)", -9999);
  output().AddVar(neutralKaonTree::lampairpchi2prev, "lampairpchi2prev", "F",
                  "Lambda pair: proton Chi2PID(2212)/npts with the residual range reversed", -9999);
  output().AddVar(neutralKaonTree::lampairpichi2pirev, "lampairpichi2pirev", "F",
                  "Lambda pair: pion Chi2PID(211)/npts with the residual range reversed", -9999);
  output().AddVar(neutralKaonTree::lampairpichi2prev, "lampairpichi2prev", "F",
                  "Lambda pair: pion Chi2PID(2212)/npts with the residual range reversed", -9999);
  output().AddVar(neutralKaonTree::lampairpnhitsdedx, "lampairpnhitsdedx", "I",
                  "Lambda pair: proton collection hits stored with a valid dE/dx", -9999);
  output().AddVar(neutralKaonTree::lampairpinhitsdedx, "lampairpinhitsdedx", "I",
                  "Lambda pair: pion collection hits stored with a valid dE/dx", -9999);
  output().AddVar(neutralKaonTree::lampairptrmdedx, "lampairptrmdedx", "F",
                  "Lambda pair: proton truncated-mean collection dE/dx [MeV/cm]", -9999);
  output().AddVar(neutralKaonTree::lampairpitrmdedx, "lampairpitrmdedx", "F",
                  "Lambda pair: pion truncated-mean collection dE/dx [MeV/cm]", -9999);
  output().AddVar(neutralKaonTree::lampairpchi2pi, "lampairpchi2pi", "F",
                  "Lambda pair: proton Chi2PID(211)/npts", -9999);
  output().AddVar(neutralKaonTree::lampairpnhitsind, "lampairpnhitsind", "I",
                  "Lambda pair: proton hits on the induction plane with more hits", -9999);
  output().AddVar(neutralKaonTree::lampairpinhitsind, "lampairpinhitsind", "I",
                  "Lambda pair: pion hits on the induction plane with more hits", -9999);
  output().AddVar(neutralKaonTree::lampairpchi2pind, "lampairpchi2pind", "F",
                  "Lambda pair: proton Chi2PID(2212)/npts on the induction plane", -9999);
  output().AddVar(neutralKaonTree::lampairpichi2pind, "lampairpichi2pind", "F",
                  "Lambda pair: pion Chi2PID(2212)/npts on the induction plane", -9999);
  output().AddVar(neutralKaonTree::lampairpichi2piind, "lampairpichi2piind", "F",
                  "Lambda pair: pion Chi2PID(211)/npts on the induction plane", -9999);
  output().AddVar(neutralKaonTree::lampairpvtxmichel, "lampairpvtxmichel", "F",
                  "Lambda pair: proton CNN Michel score around the track end", -9999);
  output().AddVar(neutralKaonTree::lampairpivtxmichel, "lampairpivtxmichel", "F",
                  "Lambda pair: pion CNN Michel score around the track end", -9999);
  output().AddVar(neutralKaonTree::lampairdoca, "lampairdoca", "F",
                  "Lambda pair: closest approach of the fitted daughter lines [cm]", -9999);
  output().AddVar(neutralKaonTree::lampairgapp, "lampairgapp", "F",
                  "Lambda pair: proton gap from the closest-approach point to its endpoint [cm]", -9999);
  output().AddVar(neutralKaonTree::lampairgappi, "lampairgappi", "F",
                  "Lambda pair: pion gap from the closest-approach point to its endpoint [cm]", -9999);
  output().AddVar(neutralKaonTree::lampaircosopen, "lampaircosopen", "F",
                  "Lambda pair: cosine of the opening angle of the fitted daughter lines", -9999);
  output().AddVar(neutralKaonTree::lampairgeompass, "lampairgeompass", "I",
                  "Lambda pair: 1 when an endpoint combination passes the geometric pairing", -9999);

  // Add neutral particle candidates variables
  const UInt_t nK0Max = neutralKaonAnalysisConstants::NMAX_K0_MICROTREE;
  neutralKaonTree::AddNeutralKaonVariables_K0Particle(output(), nK0Max);
  neutralKaonTree::AddNeutralKaonVariables_K0AlignVariants(output(), nK0Max);
  neutralKaonTree::AddNeutralKaonVariables_K0Parent(output(), nK0Max);
  neutralKaonTree::AddNeutralKaonVariables_K0CreationVtx(output(), nK0Max);
  neutralKaonTree::AddNeutralKaonVariables_K0Vtx(output(), nK0Max);
  neutralKaonTree::AddNeutralKaonVariables_K0Kinematics(output(), nK0Max);

  // Create EventDisplayData tree now that ana tree is fully defined
  // This prevents index conflicts during validation
  if (_EventDisplayDataTree) {
    _EventDisplayDataTree->InitializeTree(output());
  }

}

//********************************************************************
void neutralKaonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  pdBaseAnalysis::DefineTruthTree();
  neutralKaonTruthTree::AddNeutralKaonTruthVariables(output(), 10);
}

//********************************************************************
void neutralKaonAnalysis::FillLambdaPairDiagnostics(const AnaTrueParticleB* trueP, const AnaTrueParticleB* truePi){
//********************************************************************

  // Replays the pair requirements of pdAnnihilationUtils::CreateVerticesCommon (proton-pion mode) on the reco objects
  // of one true Lambda -> p pi- and records the first one that fails and the mask of all that fail; the pair is a
  // candidate when some neutral candidate's vertex holds reco objects of both true daughters.
  Int_t recoFlag = -1, stage = -1, failMask = -1, nHitsP = -999, nHitsPi = -999;
  Float_t dist = -999.f, pChi2P = -999.f, piChi2P = -999.f, piChi2Pi = -999.f;
  Float_t momP = -999.f, momPi = -999.f, lenP = -999.f, lenPi = -999.f;
  Int_t endProcP = -999, endProcPi = -999, recoRevP = -1, recoRevPi = -1, vtxRevP = -1, vtxRevPi = -1;
  Int_t nDedxP = -999, nDedxPi = -999;
  Float_t pChi2PRev = -999.f, piChi2PiRev = -999.f, piChi2PRev = -999.f;
  Float_t trmP = -999.f, trmPi = -999.f, pChi2Pi = -999.f;
  Int_t nHitsIndP = -999, nHitsIndPi = -999;
  Float_t pChi2PInd = -999.f, piChi2PInd = -999.f, piChi2PiInd = -999.f, michelP = -999.f, michelPi = -999.f;
  Int_t geomPass = -1;
  pdAnnihilationUtils::PairLineGeometry geometry;
  // Chi2PID per point with the residual range of the collection hits reversed (as if the track ran the other way)
  // (the residual ranges are flipped in place and restored afterwards)
  auto chi2Reversed = [](const AnaParticlePD* part, int pdg) {
    std::vector<AnaHitPD>& hits = const_cast<AnaParticlePD*>(part)->Hits[2];
    Float_t lo = 1e9f, hi = -1e9f;
    for (const AnaHitPD& hit : hits) {
      if (hit.ResidualRange < -900.f) continue;
      lo = std::min(lo, hit.ResidualRange);
      hi = std::max(hi, hit.ResidualRange);
    }
    if (hi < lo) return -999.f;
    std::vector<Float_t> saved;
    saved.reserve(hits.size());
    for (AnaHitPD& hit : hits) {
      saved.push_back(hit.ResidualRange);
      if (hit.ResidualRange > -900.f) hit.ResidualRange = hi + lo - hit.ResidualRange;
    }
    const Float_t chi2 = pdAnaUtils::Chi2PIDChi2PerHit(part, pdg);
    for (size_t i = 0; i < hits.size(); ++i) hits[i].ResidualRange = saved[i];
    return chi2;
  };
  // Truncated-mean collection dE/dx (lowest and highest 16% of the sorted valid hits dropped), as in pdAnnihilationUtils
  auto truncatedMean = [](const AnaParticlePD* part) {
    std::vector<Float_t> v;
    for (const AnaHitPD& hit : part->Hits[2]) {
      if (hit.dEdx > 0.f && hit.dEdx <= 1000.f) v.push_back(hit.dEdx);
    }
    if (v.size() < 3) return -999.f;
    std::sort(v.begin(), v.end());
    const size_t lo = static_cast<size_t>(std::lround(0.16 * v.size()));
    const size_t hi = v.size() - lo;
    double sum = 0.;
    for (size_t k = lo; k < hi; ++k) sum += v[k];
    return (hi > lo) ? static_cast<Float_t>(sum / (hi - lo)) : -999.f;
  };
  auto nValidDedx = [](const AnaParticlePD* part) {
    Int_t n = 0;
    for (const AnaHitPD& hit : part->Hits[2]) {
      if (hit.dEdx <= 1000.f && hit.dEdx != -999.f) ++n;
    }
    return n;
  };
  // 1 when the reco track starts nearer the true end than the true start
  auto recoReversed = [](const AnaParticlePD* reco, const AnaTrueParticleB* t) {
    if (reco->PositionStart[0] < -900 || reco->PositionEnd[0] < -900 || t->Position[0] < -900) return -1;
    const TVector3 ts(t->Position[0], t->Position[1], t->Position[2]);
    const TVector3 rs(reco->PositionStart[0], reco->PositionStart[1], reco->PositionStart[2]);
    const TVector3 re(reco->PositionEnd[0], reco->PositionEnd[1], reco->PositionEnd[2]);
    return ((re - ts).Mag() < (rs - ts).Mag()) ? 1 : 0;
  };
  auto trueLength = [](const AnaTrueParticleB* t) {
    const TVector3 a(t->Position[0], t->Position[1], t->Position[2]);
    const TVector3 b(t->PositionEnd[0], t->PositionEnd[1], t->PositionEnd[2]);
    return (t->Position[0] < -900 || t->PositionEnd[0] < -900) ? -999.f : static_cast<Float_t>((b - a).Mag());
  };
  if (trueP && truePi) {
    momP = trueP->Momentum;
    momPi = truePi->Momentum;
    endProcP = static_cast<Int_t>(trueP->ProcessEnd);
    endProcPi = static_cast<Int_t>(truePi->ProcessEnd);
    lenP = trueLength(trueP);
    lenPi = trueLength(truePi);
    // The reco object of each true daughter: the one with most collection hits among those matched to it
    AnaEventB& event = GetEvent();
    AnaParticlePD* recoP = nullptr;
    AnaParticlePD* recoPi = nullptr;
    for (Int_t i = 0; i < event.nParticles; ++i) {
      AnaParticlePD* reco = static_cast<AnaParticlePD*>(event.Particles[i]);
      if (!reco || !reco->TrueObject) continue;
      if (reco->TrueObject->ID == trueP->ID && (!recoP || reco->NHitsPerPlane[2] > recoP->NHitsPerPlane[2])) recoP = reco;
      if (reco->TrueObject->ID == truePi->ID && (!recoPi || reco->NHitsPerPlane[2] > recoPi->NHitsPerPlane[2])) recoPi = reco;
    }
    recoFlag = (recoP ? 0 : 1) + (recoPi ? 0 : 2);
    if (recoP && recoPi) {
      nHitsP = recoP->NHitsPerPlane[2];
      nHitsPi = recoPi->NHitsPerPlane[2];
      pChi2P = pdAnaUtils::Chi2PIDChi2PerHit(recoP, 2212);
      piChi2P = pdAnaUtils::Chi2PIDChi2PerHit(recoPi, 2212);
      piChi2Pi = pdAnaUtils::Chi2PIDChi2PerHit(recoPi, 211);
      double best = 1e9;
      const Float_t* endsP[2] = {recoP->PositionStart, recoP->PositionEnd};
      const Float_t* endsPi[2] = {recoPi->PositionStart, recoPi->PositionEnd};
      for (int ia = 0; ia < 2; ++ia) {
        for (int ib = 0; ib < 2; ++ib) {
          const Float_t* a = endsP[ia];
          const Float_t* b = endsPi[ib];
          if (a[0] < -900 || b[0] < -900) continue;
          const double d = (TVector3(a[0], a[1], a[2]) - TVector3(b[0], b[1], b[2])).Mag();
          if (d < best) {
            best = d;
            vtxRevP = ia;
            vtxRevPi = ib;
          }
        }
      }
      recoRevP = recoReversed(recoP, trueP);
      recoRevPi = recoReversed(recoPi, truePi);
      nDedxP = nValidDedx(recoP);
      nDedxPi = nValidDedx(recoPi);
      pChi2PRev = chi2Reversed(recoP, 2212);
      piChi2PiRev = chi2Reversed(recoPi, 211);
      piChi2PRev = chi2Reversed(recoPi, 2212);
      trmP = truncatedMean(recoP);
      trmPi = truncatedMean(recoPi);
      pChi2Pi = pdAnaUtils::Chi2PIDChi2PerHit(recoP, 211);
      if (recoP->NHitsPerPlane[0] > -900 || recoP->NHitsPerPlane[1] > -900)
        nHitsIndP = std::max(recoP->NHitsPerPlane[0], recoP->NHitsPerPlane[1]);
      if (recoPi->NHitsPerPlane[0] > -900 || recoPi->NHitsPerPlane[1] > -900)
        nHitsIndPi = std::max(recoPi->NHitsPerPlane[0], recoPi->NHitsPerPlane[1]);
      pChi2PInd = pdAnaUtils::Chi2PIDChi2PerHitInduction(recoP, 2212);
      piChi2PInd = pdAnaUtils::Chi2PIDChi2PerHitInduction(recoPi, 2212);
      piChi2PiInd = pdAnaUtils::Chi2PIDChi2PerHitInduction(recoPi, 211);
      michelP = recoP->vtx_CNN_michelscore;
      michelPi = recoPi->vtx_CNN_michelscore;
      if (best < 1e8) dist = static_cast<Float_t>(best);
      // Line geometry as in the geometric pairing of the vertex builder, whatever the endpoint distance
      {
        const bool combinatorics =
            !ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics") ||
            ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics") != 0;
        geomPass = pdAnnihilationUtils::BestGeometricPairing(
                       recoP, recoPi, combinatorics, ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength"),
                       ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart"),
                       pdAnnihilationUtils::ReadGeometricPairingCuts(), geometry) >= 0
                       ? 1
                       : 0;
      }

      const int minHits = ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexMinCollectionHits");
      const double maxDist = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
      const double protonBelow =
          ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDaughterProtonChi2NdfRejectBelow")
              ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDaughterProtonChi2NdfRejectBelow")
              : 100.0;
      const double pionAbove =
          ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexDaughterPionChi2NdfRejectAbove")
              ? ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexDaughterPionChi2NdfRejectAbove")
              : 25.0;
      // Hits and Chi2PID on the plane the vertex builder uses (the collection plane unless the induction fallback
      // applies, pdAnnihilationUtils::VertexDaughterHits)
      const bool fallback = ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexInductionFallback") &&
                            ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexInductionFallback") == 1;
      const Int_t effHitsP = pdAnnihilationUtils::VertexDaughterHits(recoP, fallback, minHits);
      const Int_t effHitsPi = pdAnnihilationUtils::VertexDaughterHits(recoPi, fallback, minHits);
      const Float_t effPChi2P = pdAnnihilationUtils::VertexDaughterChi2PerHit(recoP, 2212, fallback, minHits);
      const Float_t effPiChi2P = pdAnnihilationUtils::VertexDaughterChi2PerHit(recoPi, 2212, fallback, minHits);
      const Float_t effPiChi2Pi = pdAnnihilationUtils::VertexDaughterChi2PerHit(recoPi, 211, fallback, minHits);
      const bool pIsProtonLike = effPChi2P > -900.f && effPChi2P < protonBelow;
      const bool piIsProtonLike = effPiChi2P > -900.f && effPiChi2P < protonBelow;
      // Pion veto as in pdAnnihilationUtils (IsProtonLikeAndNotPionLike is file-local there)
      const bool piPassesVeto =
          !(effPiChi2P > -900.f && (effPiChi2P < protonBelow || (effPiChi2Pi > -900.f && effPiChi2Pi > pionAbove)));

      // A candidate by truth (either reco object of each daughter), as in the analysis signal definition
      bool isCandidate = false;
      const ToyBoxNeutralKaon& nkBox = static_cast<const ToyBoxNeutralKaon&>(box());
      for (const AnaNeutralParticlePD* cand : nkBox.neutralParticleCandidates) {
        const AnaAnnihilationVertexPD* v = cand ? cand->AnnihilationVertex : nullptr;
        if (!v || v->Particles.size() < 2 || !v->Particles[0] || !v->Particles[1]) continue;
        const auto* t0 = v->Particles[0]->TrueObject;
        const auto* t1 = v->Particles[1]->TrueObject;
        if (!t0 || !t1) continue;
        if ((t0->ID == trueP->ID && t1->ID == truePi->ID) || (t0->ID == truePi->ID && t1->ID == trueP->ID)) {
          isCandidate = true;
          break;
        }
      }
      const bool failHits = recoP->PositionStart[0] < -900 || recoPi->PositionStart[0] < -900 || effHitsP <= minHits ||
                            effHitsPi <= minHits;
      const bool failParent = recoP->ParentID != recoPi->ParentID;
      const bool failPID = !pIsProtonLike || piIsProtonLike || !piPassesVeto;
      const bool failDist = !(dist > -900.f) || dist > maxDist;
      failMask = (failHits ? 1 : 0) | (failParent ? 2 : 0) | (failPID ? 4 : 0) | (failDist ? 8 : 0);
      if (isCandidate) stage = 0;
      else if (failHits) stage = 1;
      else if (failParent) stage = 2;
      else if (failPID) stage = 3;
      else if (failDist) stage = 4;
      else stage = 5;
    }
  }
  output().FillVar(neutralKaonTree::lampairrecoflag, recoFlag);
  output().FillVar(neutralKaonTree::lampairstage, stage);
  output().FillVar(neutralKaonTree::lampairfailmask, failMask);
  output().FillVar(neutralKaonTree::lampairdist, dist);
  output().FillVar(neutralKaonTree::lampairnhitsp, nHitsP);
  output().FillVar(neutralKaonTree::lampairnhitspi, nHitsPi);
  output().FillVar(neutralKaonTree::lampairpchi2p, pChi2P);
  output().FillVar(neutralKaonTree::lampairpichi2p, piChi2P);
  output().FillVar(neutralKaonTree::lampairpichi2pi, piChi2Pi);
  output().FillVar(neutralKaonTree::lampairtruemomp, momP);
  output().FillVar(neutralKaonTree::lampairtruemompi, momPi);
  output().FillVar(neutralKaonTree::lampairtruelengthp, lenP);
  output().FillVar(neutralKaonTree::lampairtruelengthpi, lenPi);
  output().FillVar(neutralKaonTree::lampairptrueendproc, endProcP);
  output().FillVar(neutralKaonTree::lampairpitrueendproc, endProcPi);
  output().FillVar(neutralKaonTree::lampairprecoreversed, recoRevP);
  output().FillVar(neutralKaonTree::lampairpirecoreversed, recoRevPi);
  output().FillVar(neutralKaonTree::lampairpvtxreversed, vtxRevP);
  output().FillVar(neutralKaonTree::lampairpivtxreversed, vtxRevPi);
  output().FillVar(neutralKaonTree::lampairpchi2prev, pChi2PRev);
  output().FillVar(neutralKaonTree::lampairpichi2pirev, piChi2PiRev);
  output().FillVar(neutralKaonTree::lampairpichi2prev, piChi2PRev);
  output().FillVar(neutralKaonTree::lampairpnhitsdedx, nDedxP);
  output().FillVar(neutralKaonTree::lampairpinhitsdedx, nDedxPi);
  output().FillVar(neutralKaonTree::lampairptrmdedx, trmP);
  output().FillVar(neutralKaonTree::lampairpitrmdedx, trmPi);
  output().FillVar(neutralKaonTree::lampairpchi2pi, pChi2Pi);
  output().FillVar(neutralKaonTree::lampairpnhitsind, nHitsIndP);
  output().FillVar(neutralKaonTree::lampairpinhitsind, nHitsIndPi);
  output().FillVar(neutralKaonTree::lampairpchi2pind, pChi2PInd);
  output().FillVar(neutralKaonTree::lampairpichi2pind, piChi2PInd);
  output().FillVar(neutralKaonTree::lampairpichi2piind, piChi2PiInd);
  output().FillVar(neutralKaonTree::lampairpvtxmichel, michelP);
  output().FillVar(neutralKaonTree::lampairpivtxmichel, michelPi);
  output().FillVar(neutralKaonTree::lampairdoca, geometry.doca);
  output().FillVar(neutralKaonTree::lampairgapp, geometry.gap1);
  output().FillVar(neutralKaonTree::lampairgappi, geometry.gap2);
  output().FillVar(neutralKaonTree::lampaircosopen, geometry.cosOpening);
  output().FillVar(neutralKaonTree::lampairgeompass, geomPass);
}

//********************************************************************
void neutralKaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  if (addBase) pdBaseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis (only once)
  standardPDTree::FillStandardVariables_EventInfo(output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(output(), GetSpill().Beam);

  // Lambda -> p pi- census. The minitrees often drop the Lambda itself but keep reco-matched particles, so a
  // Lambda decay is recognised from its daughters: a true p and pi- with the same parent, both created by a decay.
  {
    std::set<Int_t> recoMatchedIDs;
    for (size_t i = 0; i < GetBunch().Particles.size(); ++i) {
      const AnaParticleB* reco = GetBunch().Particles[i];
      if (reco && reco->TrueObject) recoMatchedIDs.insert(reco->TrueObject->ID);
    }
    std::map<Int_t, std::pair<const AnaTrueParticleB*, const AnaTrueParticleB*> > byParent;  // (p, pi-)
    for (size_t i = 0; i < GetSpill().TrueParticles.size(); ++i) {
      const AnaTrueParticleB* t = GetSpill().TrueParticles[i];
      if (!t || t->ProcessStart != AnaTrueParticleB::Decay) continue;
      if (t->PDG == 2212) byParent[t->ParentID].first = t;
      if (t->PDG == -211) byParent[t->ParentID].second = t;
    }
    Int_t nPairs = 0, nPairsReco = 0;
    const AnaTrueParticleB* diagP = nullptr;
    const AnaTrueParticleB* diagPi = nullptr;
    for (const auto& kv : byParent) {
      if (!kv.second.first || !kv.second.second) continue;
      ++nPairs;
      const bool bothReco = recoMatchedIDs.count(kv.second.first->ID) && recoMatchedIDs.count(kv.second.second->ID);
      if (bothReco) ++nPairsReco;
      if (!diagP || (bothReco && nPairsReco == 1)) {
        diagP = kv.second.first;
        diagPi = kv.second.second;
      }
    }
    output().FillVar(neutralKaonTree::ntruelambdapairs, nPairs);
    output().FillVar(neutralKaonTree::ntruelambdapairsreco, nPairsReco);
    FillLambdaPairDiagnostics(diagP, diagPi);
  }

  // Get neutral kaon box for accessing candidates
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());

  // Fill individual candidate data
  if(neutralKaonBox.neutralParticleCandidates.size() > 0){
    for(size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++){
      neutralKaonTree::FillNeutralKaonVariables(output(), neutralKaonBox.neutralParticleCandidates[i], GetEvent(),
                                                neutralKaonBox.nAnnihilationVerticesBeforeFiltering,
                                                neutralKaonBox.nAnnihilationVerticesAfterFiltering,
                                                GetSpill().Beam, i);
      output().IncrementCounter(neutralKaonTree::nk0);
    }
  }

  // Save event display data for all events
  if (_EventDisplayDataTree) {
    _EventDisplayDataTree->FillEventDisplayData(output(), GetEvent(), neutralKaonBox);
  }

//  // Save event display data only for events with a signal neutral candidate (signal==1)
//  if (_EventDisplayDataTree && !neutralKaonBox.neutralParticleCandidates.empty()) {
//    bool hasSignalCandidate = false;

//    if (anaUtils::_categ && anaUtils::_categ->HasCategory("signal")) {
//      TrackCategoryDefinition& signalCategory = anaUtils::_categ->GetCategory("signal");
//      for (size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); ++i) {
//        if (signalCategory.GetObjectCode(1, static_cast<Int_t>(i)) != -999) {
//          hasSignalCandidate = true;
//          break;
//        }
//      }
//    }

//    if (hasSignalCandidate) {
//      _EventDisplayDataTree->FillEventDisplayData(output(), GetEvent(), neutralKaonBox);
//    }
//  }
}

//********************************************************************
void neutralKaonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

   // Fill the common variables
  if (addBase) pdBaseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
bool neutralKaonAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  (void) vtx; // to avoid warning for unused vtx variable

  // fill it allways for the moment
  return true;
}
//********************************************************************
bool neutralKaonAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* part){
//********************************************************************
  if (!part) return false;
  // PDG of the true neutral the truth tree is filled for (default K0S; e.g. 3122 to study Lambda -> p pi-)
  static const Int_t truthTreePDG = ND::params().HasParameter("neutralKaonAnalysis.TruthTreePDG")
                                        ? ND::params().GetParameterI("neutralKaonAnalysis.TruthTreePDG")
                                        : 310;
  if(part->PDG != truthTreePDG) return false;
  return true;
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************
  pdBaseAnalysis::FillTruthTree(part);

  AnaTrueParticlePD* daughter1True = nullptr;
  AnaTrueParticlePD* daughter2True = nullptr;
  AnaTrueParticlePD* parentTrue = nullptr;

  for(int i = 0; i < GetSpill().TrueParticles.size(); ++i){
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[i]);
    if(!truePart) continue;
    if(part.Daughters.size() > 0 && truePart->ID == part.Daughters[0]) daughter1True = truePart;
    if(part.Daughters.size() > 1 && truePart->ID == part.Daughters[1]) daughter2True = truePart;
    if(truePart->ID == part.ParentID) parentTrue = truePart;
  }

  auto findRecoFromTrue = [&](const AnaTrueParticlePD* truePart, bool includeBeam) -> AnaParticlePD* {
    if(!truePart) return nullptr;

    for(size_t i = 0; i < GetBunch().Particles.size(); ++i){
      AnaParticlePD* recoPart = static_cast<AnaParticlePD*>(GetBunch().Particles[i]);
      if(!recoPart || !recoPart->TrueObject) continue;
      AnaTrueParticlePD* recoTruePart = static_cast<AnaTrueParticlePD*>(recoPart->TrueObject);
      if(recoTruePart && recoTruePart->ID == truePart->ID) return recoPart;
    }

    if(includeBeam){
      AnaParticlePD* beamParticle = static_cast<AnaParticlePD*>(static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle);
      if(beamParticle && beamParticle->TrueObject){
        AnaTrueParticlePD* trueBeam = static_cast<AnaTrueParticlePD*>(beamParticle->TrueObject);
        if(trueBeam && trueBeam->ID == truePart->ID) return beamParticle;
      }
    }

    return nullptr;
  };

  AnaParticlePD* daughter1Reco = findRecoFromTrue(daughter1True, false);
  AnaParticlePD* daughter2Reco = findRecoFromTrue(daughter2True, false);

  auto trueSubDaughterRecoStats = [&](const AnaTrueParticlePD* mother) -> std::pair<Int_t, Int_t> {
    if (!mother) return {-999, -999};
    Int_t nQual = 0;
    Int_t hitsTot = 0;
    for (Int_t childId : mother->Daughters) {
      AnaTrueParticlePD* childTrue = nullptr;
      for (int ti = 0; ti < GetSpill().TrueParticles.size(); ++ti) {
        AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[ti]);
        if (tp && tp->ID == childId) {
          childTrue = tp;
          break;
        }
      }
      if (!childTrue) continue;
      AnaParticlePD* childReco = findRecoFromTrue(childTrue, false);
      if (!HasRecoObjectForTruthFlags(childReco)) continue;
      ++nQual;
      const Int_t nh = RecoNHitsOrSentinel(childReco);
      if (nh >= 0) hitsTot += nh;
    }
    return {nQual, hitsTot};
  };

  auto trueSubDaughterRecoTrueEnergySum = [&](const AnaTrueParticlePD* mother) -> Float_t {
    if(!mother) return static_cast<Float_t>(-999.0);
    double sumGeV = 0.0;
    for(Int_t childId : mother->Daughters){
      AnaTrueParticlePD* childTrue = nullptr;
      for(int ti = 0; ti < GetSpill().TrueParticles.size(); ++ti){
        AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[ti]);
        if(tp && tp->ID == childId){
          childTrue = tp;
          break;
        }
      }
      if(!childTrue) continue;
      AnaParticlePD* childReco = findRecoFromTrue(childTrue, false);
      if(!HasRecoObjectForTruthFlags(childReco)) continue;
      const Float_t eGeV = TrueTotalEnergyGeV(childTrue, false);
      if(std::isfinite(eGeV) && eGeV >= 0.f && eGeV < 1.e4f) sumGeV += static_cast<double>(eGeV);
    }
    return static_cast<Float_t>(sumGeV);
  };

  const std::pair<Int_t, Int_t> daughter1TrueSubReco = trueSubDaughterRecoStats(daughter1True);
  const std::pair<Int_t, Int_t> daughter2TrueSubReco = trueSubDaughterRecoStats(daughter2True);

  AnaParticlePD* parentReco = findRecoFromTrue(parentTrue, true);

  const bool parentHasReco = HasRecoObjectForTruthFlags(parentReco);
  const bool daughter1HasReco = HasRecoObjectForTruthFlags(daughter1Reco);
  const bool daughter2HasReco = HasRecoObjectForTruthFlags(daughter2Reco);
  const Float_t daughterStartDistReco = RecoDaughtersStartDistance(daughter1Reco, daughter2Reco);
  const Int_t daughterRecoParentIdDifference =
    (daughter1Reco && daughter2Reco) ? (daughter1Reco->ParentID - daughter2Reco->ParentID) : -999;
  const Int_t daughter1RecoParentVsTrueParentRecoIdDiff =
    RecoDaughterParentUidMinusTrueParentRecoUid(GetBunch(), daughter1Reco, parentReco);
  const Int_t daughter2RecoParentVsTrueParentRecoIdDiff =
    RecoDaughterParentUidMinusTrueParentRecoUid(GetBunch(), daughter2Reco, parentReco);
  const Float_t parentRecoTrueStartDistance = RecoVsTrueStartDistance(parentReco, parentTrue);
  const Float_t parentRecoTrueEndDistance = RecoVsTrueEndDistance(parentReco, parentTrue);
  const Float_t daughter1RecoTrueStartDistance = RecoVsTrueStartDistance(daughter1Reco, daughter1True);
  const Float_t daughter1RecoTrueEndDistance = RecoVsTrueEndDistance(daughter1Reco, daughter1True);
  const Float_t daughter2RecoTrueStartDistance = RecoVsTrueStartDistance(daughter2Reco, daughter2True);
  const Float_t daughter2RecoTrueEndDistance = RecoVsTrueEndDistance(daughter2Reco, daughter2True);

  const double annihilationVertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");

  Float_t daughterPairFitMinDistance = static_cast<Float_t>(-999.0);
  if(daughter1Reco && daughter2Reco){
    const Float_t dStart = RecoDaughtersStartDistance(daughter1Reco, daughter2Reco);
    if(IsValidDistanceValue(dStart) && dStart <= static_cast<Float_t>(annihilationVertexRadius)){
      daughterPairFitMinDistance = RecoPairFitMinimumDistance(daughter1Reco, daughter2Reco,
                                                              trackFitLength, trackFitDistanceFromStart);
    }
  }

  Float_t daughter1OtherPairFitMinDistance = static_cast<Float_t>(-999.0);
  Float_t daughter2OtherPairFitMinDistance = static_cast<Float_t>(-999.0);
  Float_t otherPairsFitMinDistanceGlobal = static_cast<Float_t>(-999.0);
  AnaParticlePD* daughter1OtherArgminPartner = nullptr;

  auto updateMin = [](Float_t& currentMin, Float_t candidate){
    if(!IsValidDistanceValue(candidate)) return;
    if(!IsValidDistanceValue(currentMin) || candidate < currentMin) currentMin = candidate;
  };

  for(size_t i = 0; i < GetBunch().Particles.size(); ++i){
    AnaParticlePD* other = static_cast<AnaParticlePD*>(GetBunch().Particles[i]);
    if(!other) continue;

    if(daughter1Reco && other != daughter1Reco && other != daughter2Reco){
      const Float_t dStart = RecoDaughtersStartDistance(daughter1Reco, other);
      if(IsValidDistanceValue(dStart) && dStart <= static_cast<Float_t>(annihilationVertexRadius)){
        const Float_t dMin = RecoPairFitMinimumDistance(daughter1Reco, other,
                                                        trackFitLength, trackFitDistanceFromStart);
        if(IsValidDistanceValue(dMin)){
          if(!IsValidDistanceValue(daughter1OtherPairFitMinDistance) || dMin < daughter1OtherPairFitMinDistance){
            daughter1OtherPairFitMinDistance = dMin;
            daughter1OtherArgminPartner = other;
          }
        }
      }
    }

    if(daughter2Reco && other != daughter2Reco && other != daughter1Reco){
      const Float_t dStart = RecoDaughtersStartDistance(daughter2Reco, other);
      if(IsValidDistanceValue(dStart) && dStart <= static_cast<Float_t>(annihilationVertexRadius)){
        const Float_t dMin = RecoPairFitMinimumDistance(daughter2Reco, other,
                                                        trackFitLength, trackFitDistanceFromStart);
        updateMin(daughter2OtherPairFitMinDistance, dMin);
      }
    }
  }

  updateMin(otherPairsFitMinDistanceGlobal, daughter1OtherPairFitMinDistance);
  updateMin(otherPairsFitMinDistanceGlobal, daughter2OtherPairFitMinDistance);

  const Int_t daughter1OtherMinSepTruePdg = TruePdgOfReco(daughter1OtherArgminPartner);
  const Int_t vertexMinIsK0DaughterPair =
    VertexGlobalMinIsDauPairTag(daughterPairFitMinDistance, daughter1OtherPairFitMinDistance,
                                daughter2OtherPairFitMinDistance);

  const Float_t daughter1RecoProtonChi2Ndf = RecoProtonChi2Ndf(daughter1Reco);
  const Float_t daughter1RecoKaonChi2Ndf = RecoKaonChi2Ndf(daughter1Reco);
  const Float_t daughter1RecoPionChi2Ndf = RecoPionChi2Ndf(daughter1Reco);
  const Float_t daughter2RecoProtonChi2Ndf = RecoProtonChi2Ndf(daughter2Reco);
  const Float_t daughter2RecoKaonChi2Ndf = RecoKaonChi2Ndf(daughter2Reco);
  const Float_t daughter2RecoPionChi2Ndf = RecoPionChi2Ndf(daughter2Reco);

  const Float_t k0TrueMomentum = TrueScalarMomentum(&part, false);
  const Float_t parentTrueStartMomentum = TrueScalarMomentum(parentTrue, false);
  const Float_t parentTrueEndMomentum = TrueScalarMomentum(parentTrue, true);
  const Float_t daughter1TrueStartMomentum = TrueScalarMomentum(daughter1True, false);
  const Float_t daughter1TrueEndMomentum = TrueScalarMomentum(daughter1True, true);
  const Float_t daughter2TrueStartMomentum = TrueScalarMomentum(daughter2True, false);
  const Float_t daughter2TrueEndMomentum = TrueScalarMomentum(daughter2True, true);

  const Float_t k0TrueEnergy = TrueTotalEnergyGeV(&part, false);
  const Float_t daughter1TrueEnergy = TrueTotalEnergyGeV(daughter1True, false);
  const Float_t daughter2TrueEnergy = TrueTotalEnergyGeV(daughter2True, false);
  const Float_t daughter1TrueSubRecoEnergySum = trueSubDaughterRecoTrueEnergySum(daughter1True);
  const Float_t daughter2TrueSubRecoEnergySum = trueSubDaughterRecoTrueEnergySum(daughter2True);
  const Int_t isK0Decay = (part.ProcessEnd == AnaTrueParticleB::Decay) ? 1 : 0;
  const Int_t k0NTrueDaughters = static_cast<Int_t>(part.Daughters.size());

  const Int_t k0TruePdg = part.PDG;
  const Int_t parentTruePdg = parentTrue ? parentTrue->PDG : -999;
  const Int_t daughter1TruePdg = daughter1True ? daughter1True->PDG : -999;
  const Int_t daughter2TruePdg = daughter2True ? daughter2True->PDG : -999;
  const Int_t isK0Charged =
    ((daughter1TruePdg == 211 && daughter2TruePdg == -211) ||
     (daughter1TruePdg == -211 && daughter2TruePdg == 211)) ? 1 : 0;
  const Int_t isK0Neutral =
    ((daughter1TruePdg == 111 && daughter2TruePdg == 111)) ? 1 : 0;
  auto isPi0ExactlyToTwoGamma = [&](const AnaTrueParticlePD* pi0) -> Int_t {
    if(!pi0 || pi0->PDG != 111) return 0;
    if(pi0->Daughters.size() != 2) return 0;

    AnaTrueParticlePD* gamma1 = nullptr;
    AnaTrueParticlePD* gamma2 = nullptr;
    for(int ti = 0; ti < GetSpill().TrueParticles.size(); ++ti){
      AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[ti]);
      if(!tp) continue;
      if(tp->ID == pi0->Daughters[0]) gamma1 = tp;
      if(tp->ID == pi0->Daughters[1]) gamma2 = tp;
    }

    if(!gamma1 || !gamma2) return 0;
    return (gamma1->PDG == 22 && gamma2->PDG == 22) ? 1 : 0;
  };
  const Int_t k0Pi01TwoGamma = isPi0ExactlyToTwoGamma(daughter1True);
  const Int_t k0Pi02TwoGamma = isPi0ExactlyToTwoGamma(daughter2True);
  auto gammaHasValidPandoraRecoObject = [&](const AnaTrueParticlePD* pi0, Int_t gammaDaughterIndex) -> Int_t {
    if(!pi0 || pi0->PDG != 111) return 0;
    if(gammaDaughterIndex < 0 || gammaDaughterIndex >= static_cast<Int_t>(pi0->Daughters.size())) return 0;
    const Int_t gammaId = pi0->Daughters[gammaDaughterIndex];

    AnaTrueParticlePD* gammaTrue = nullptr;
    for(int ti = 0; ti < GetSpill().TrueParticles.size(); ++ti){
      AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[ti]);
      if(tp && tp->ID == gammaId){
        gammaTrue = tp;
        break;
      }
    }
    if(!gammaTrue || gammaTrue->PDG != 22) return 0;

    AnaParticlePD* gammaReco = findRecoFromTrue(gammaTrue, false);
    return (gammaReco && gammaReco->isPandora && HasRecoObjectForTruthFlags(gammaReco)) ? 1 : 0;
  };
  const Int_t k0Dau1Gamma1HasRecoObject = gammaHasValidPandoraRecoObject(daughter1True, 0);
  const Int_t k0Dau1Gamma2HasRecoObject = gammaHasValidPandoraRecoObject(daughter1True, 1);
  const Int_t k0Dau2Gamma1HasRecoObject = gammaHasValidPandoraRecoObject(daughter2True, 0);
  const Int_t k0Dau2Gamma2HasRecoObject = gammaHasValidPandoraRecoObject(daughter2True, 1);

  const Float_t daughter1TrueLength = TrueLengthVertexSeparationCm(daughter1True);
  const Float_t daughter1RecoLength = RecoLengthOrSentinel(daughter1Reco);
  const Float_t daughter2TrueLength = TrueLengthVertexSeparationCm(daughter2True);
  const Float_t daughter2RecoLength = RecoLengthOrSentinel(daughter2Reco);
  const Float_t k0TrueLength = TrueLengthVertexSeparationCm(&part);
  const Int_t daughter1RecoNHits = RecoNHitsOrSentinel(daughter1Reco);
  const Int_t daughter2RecoNHits = RecoNHitsOrSentinel(daughter2Reco);
  const Int_t parentRecoNHits = RecoNHitsOrSentinel(parentReco);
  const Float_t parentTrueLength = TrueLengthVertexSeparationCm(parentTrue);
  const Float_t parentRecoLength = RecoLengthOrSentinel(parentReco);

  const Int_t parentRecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(parentReco);
  const Int_t daughter1RecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(daughter1Reco);
  const Int_t daughter2RecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(daughter2Reco);

  const Int_t parentIsPandoraBeam =
    (parentReco && parentReco->isPandora) ? 1 : 0;

  // Among true K0 parent's direct true daughters: exactly one has valid reco, and that particle is a proton (2212).
  Int_t parentExactlyOneRecoProtonNearK0Creation = 0;
  if(parentTrue){
    int nSiblingsWithReco = 0;
    AnaTrueParticlePD* soleSiblingTrueWithReco = nullptr;
    for(Int_t childId : parentTrue->Daughters){
      AnaTrueParticlePD* childTrue = nullptr;
      for(int ti = 0; ti < GetSpill().TrueParticles.size(); ++ti){
        AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[ti]);
        if(tp && tp->ID == childId){
          childTrue = tp;
          break;
        }
      }
      if(!childTrue) continue;
      AnaParticlePD* siblingReco = findRecoFromTrue(childTrue, false);
      if(!HasRecoObjectForTruthFlags(siblingReco)) continue;
      ++nSiblingsWithReco;
      soleSiblingTrueWithReco = childTrue;
    }
    if(nSiblingsWithReco == 1 && soleSiblingTrueWithReco && soleSiblingTrueWithReco->PDG == 2212)
      parentExactlyOneRecoProtonNearK0Creation = 1;
  }

  neutralKaonTruthTree::FillNeutralKaonTruthVariables(output(),
                                                      parentHasReco,
                                                      daughter1HasReco,
                                                      daughter2HasReco,
                                                      daughterStartDistReco,
                                                      daughterRecoParentIdDifference,
                                                      daughter1RecoParentVsTrueParentRecoIdDiff,
                                                      daughter2RecoParentVsTrueParentRecoIdDiff,
                                                      parentRecoTrueStartDistance,
                                                      parentRecoTrueEndDistance,
                                                      daughter1RecoTrueStartDistance,
                                                      daughter1RecoTrueEndDistance,
                                                      daughter2RecoTrueStartDistance,
                                                      daughter2RecoTrueEndDistance,
                                                      daughterPairFitMinDistance,
                                                      daughter1OtherPairFitMinDistance,
                                                      daughter2OtherPairFitMinDistance,
                                                      otherPairsFitMinDistanceGlobal,
                                                      daughter1RecoProtonChi2Ndf,
                                                      daughter1RecoKaonChi2Ndf,
                                                      daughter1RecoPionChi2Ndf,
                                                      daughter2RecoProtonChi2Ndf,
                                                      daughter2RecoKaonChi2Ndf,
                                                      daughter2RecoPionChi2Ndf,
                                                      k0TrueMomentum,
                                                      parentTrueStartMomentum,
                                                      parentTrueEndMomentum,
                                                      daughter1TrueStartMomentum,
                                                      daughter1TrueEndMomentum,
                                                      daughter2TrueStartMomentum,
                                                      daughter2TrueEndMomentum,
                                                      k0TrueEnergy,
                                                      daughter1TrueEnergy,
                                                      daughter2TrueEnergy,
                                                      daughter1TrueSubRecoEnergySum,
                                                      daughter2TrueSubRecoEnergySum,
                                                      isK0Decay,
                                                      k0NTrueDaughters,
                                                      isK0Charged,
                                                      isK0Neutral,
                                                      k0Pi01TwoGamma,
                                                      k0Pi02TwoGamma,
                                                      k0Dau1Gamma1HasRecoObject,
                                                      k0Dau1Gamma2HasRecoObject,
                                                      k0Dau2Gamma1HasRecoObject,
                                                      k0Dau2Gamma2HasRecoObject,
                                                      k0TruePdg,
                                                      parentTruePdg,
                                                      daughter1TruePdg,
                                                      daughter2TruePdg,
                                                      daughter1TrueLength,
                                                      daughter1RecoLength,
                                                      daughter2TrueLength,
                                                      daughter2RecoLength,
                                                      k0TrueLength,
                                                      daughter1RecoNHits,
                                                      daughter2RecoNHits,
                                                      daughter1TrueSubReco.first,
                                                      daughter1TrueSubReco.second,
                                                      daughter2TrueSubReco.first,
                                                      daughter2TrueSubReco.second,
                                                      parentRecoNHits,
                                                      parentTrueLength,
                                                      parentRecoLength,
                                                      daughter1OtherMinSepTruePdg,
                                                      vertexMinIsK0DaughterPair,
                                                      parentRecoZStartGreaterThanEnd,
                                                      daughter1RecoZStartGreaterThanEnd,
                                                      daughter2RecoZStartGreaterThanEnd,
                                                      parentIsPandoraBeam,
                                                      parentExactlyOneRecoProtonNearK0Creation);
  // Line geometry of the reco daughters as in the geometric pairing, whatever their start distance
  {
    pdAnnihilationUtils::PairLineGeometry geometry;
    Int_t geomPass = -1;
    if(daughter1Reco && daughter2Reco){
      const bool combinatorics =
        !ND::params().HasParameter("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics") ||
        ND::params().GetParameterI("neutralKaonAnalysis.AnnihilationVertexUseEndpointCombinatorics") != 0;
      geomPass = (pdAnnihilationUtils::BestGeometricPairing(daughter1Reco, daughter2Reco, combinatorics,
                                                            trackFitLength, trackFitDistanceFromStart,
                                                            pdAnnihilationUtils::ReadGeometricPairingCuts(),
                                                            geometry) >= 0) ? 1 : 0;
    }
    neutralKaonTruthTree::FillNeutralKaonTruthPairGeometry(output(), geometry.doca, geometry.gap1, geometry.gap2,
                                                           geometry.cosOpening, geomPass);
  }
  output().FillVectorVar(neutralKaonTruthTree::k0trueendx, part.PositionEnd[0]);
  output().FillVectorVar(neutralKaonTruthTree::k0trueendy, part.PositionEnd[1]);
  output().FillVectorVar(neutralKaonTruthTree::k0trueendz, part.PositionEnd[2]);
  output().IncrementCounter(neutralKaonTruthTree::ntruek0);
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables
  pdBaseAnalysis::FillTruthTreeBase(vtx);
}

// bool neutralKaonAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* part){
//   return true;
// }

// void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//   // Fill the common variables
//   pdBaseAnalysis::FillTruthTree();
// }

//********************************************************************
void neutralKaonAnalysis::FillCategories(){
//********************************************************************

  // For the candidate
  if(box().MainTrack)
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,"");

  // For the beam track
  AnaParticleB* beam = static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle;
  if(beam)anaUtils::FillCategories(&GetEvent(), beam,"beam");

  // For neutral particle candidates
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());
  if(neutralKaonBox.neutralParticleCandidates.size() > 0){
    for(size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++){
      neutralKaonAnaUtils::FillSignalCandidateCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonAnaUtils::FillVertexCandidateCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
    }
  }

}