#include "neutralKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "neutralKaonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"
#include "neutralKaonTruthTree.hxx"
#include "ToyBoxNeutralKaon.hxx"

#include "pdAnalysisUtils.hxx"
#include "standardPDTree.hxx"
#include <sstream>

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

  Float_t TrueLengthOrSentinel(const AnaTrueParticlePD* truePart){
    if(!truePart) return static_cast<Float_t>(-999.0);
    const Float_t L = truePart->Length;
    if(!std::isfinite(L) || L < 0.f) return static_cast<Float_t>(-999.0);
    return L;
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

  Float_t RecoChi2NdfFromPidPlane(const AnaParticlePD* recoPart, int plane, int pidIndex){
    if(!recoPart || plane < 0 || plane > 2) return static_cast<Float_t>(-999.0);
    const Float_t ndf = recoPart->PID[plane][0];
    const Float_t chi2 = recoPart->PID[plane][pidIndex];
    if(ndf <= 0.f || !std::isfinite(ndf) || !std::isfinite(chi2) || chi2 < 0.f) return static_cast<Float_t>(-999.0);
    return chi2 / ndf;
  }

  Float_t RecoChi2NdfFromPidBestPlane(const AnaParticlePD* recoPart, int pidIndex){
    if(!recoPart) return static_cast<Float_t>(-999.0);
    for(int plane = 2; plane >= 0; --plane){
      const Float_t v = RecoChi2NdfFromPidPlane(recoPart, plane, pidIndex);
      if(std::isfinite(v) && v > -900.f) return v;
    }
    return static_cast<Float_t>(-999.0);
  }

  Float_t RecoProtonChi2Ndf(const AnaParticlePD* recoPart){
    if(!recoPart) return static_cast<Float_t>(-999.0);
    if(recoPart->Chi2ndf > 0.f && recoPart->Chi2Proton > 0.f &&
       std::isfinite(recoPart->Chi2ndf) && std::isfinite(recoPart->Chi2Proton))
      return recoPart->Chi2Proton / recoPart->Chi2ndf;
    return RecoChi2NdfFromPidBestPlane(recoPart, 3);
  }

  Float_t RecoKaonChi2Ndf(const AnaParticlePD* recoPart){
    return RecoChi2NdfFromPidBestPlane(recoPart, 4);
  }

  Float_t RecoPionChi2Ndf(const AnaParticlePD* recoPart){
    return RecoChi2NdfFromPidBestPlane(recoPart, 5);
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
void neutralKaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  if (addBase) pdBaseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis (only once)
  standardPDTree::FillStandardVariables_EventInfo(output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(output(), GetSpill().Beam);

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
  if(part->PDG != 310) return false;
  if(part->ProcessEnd != AnaTrueParticleB::Decay) return false;
  if(part->Daughters.size() != 2) return false;

  AnaTrueParticlePD* daughter1 = nullptr;
  AnaTrueParticlePD* daughter2 = nullptr;

  for(int i = 0; i < GetSpill().TrueParticles.size(); ++i){
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[i]);
    if(!truePart) continue;
    if(truePart->ID == part->Daughters[0]) daughter1 = truePart;
    if(truePart->ID == part->Daughters[1]) daughter2 = truePart;
  }

  if(!daughter1 || !daughter2) return false;

  const bool isPiPlusPiMinus =
    ((daughter1->PDG == 211 && daughter2->PDG == -211) ||
     (daughter1->PDG == -211 && daughter2->PDG == 211));

  return isPiPlusPiMinus;
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

  const Int_t k0TruePdg = part.PDG;
  const Int_t parentTruePdg = parentTrue ? parentTrue->PDG : -999;
  const Int_t daughter1TruePdg = daughter1True ? daughter1True->PDG : -999;
  const Int_t daughter2TruePdg = daughter2True ? daughter2True->PDG : -999;

  const Float_t daughter1TrueLength = TrueLengthOrSentinel(daughter1True);
  const Float_t daughter1RecoLength = RecoLengthOrSentinel(daughter1Reco);
  const Float_t daughter2TrueLength = TrueLengthOrSentinel(daughter2True);
  const Float_t daughter2RecoLength = RecoLengthOrSentinel(daughter2Reco);
  const Float_t k0TrueLength = TrueLengthOrSentinel(&part);
  const Int_t daughter1RecoNHits = RecoNHitsOrSentinel(daughter1Reco);
  const Int_t daughter2RecoNHits = RecoNHitsOrSentinel(daughter2Reco);
  const Int_t parentRecoNHits = RecoNHitsOrSentinel(parentReco);
  const Float_t parentTrueLength = TrueLengthOrSentinel(parentTrue);
  const Float_t parentRecoLength = RecoLengthOrSentinel(parentReco);

  const Int_t parentRecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(parentReco);
  const Int_t daughter1RecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(daughter1Reco);
  const Int_t daughter2RecoZStartGreaterThanEnd = RecoZStartGreaterThanEndFlag(daughter2Reco);

  const Int_t parentIsPandoraBeam =
    (parentReco && parentReco->isPandora) ? 1 : 0;

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
                                                      parentRecoNHits,
                                                      parentTrueLength,
                                                      parentRecoLength,
                                                      daughter1OtherMinSepTruePdg,
                                                      vertexMinIsK0DaughterPair,
                                                      parentRecoZStartGreaterThanEnd,
                                                      daughter1RecoZStartGreaterThanEnd,
                                                      daughter2RecoZStartGreaterThanEnd,
                                                      parentIsPandoraBeam);
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