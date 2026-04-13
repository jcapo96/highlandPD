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
#include <iostream>
#include <unordered_map>

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

#include "CategoryManager.hxx"

#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"

#include "baseToyMaker.hxx"

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
  neutralKaonTruthTree::AddNeutralKaonParentTruthVariables(output(), 10);
  neutralKaonTruthTree::AddNeutralKaonDaughter1TruthVariables(output(), 10);
  neutralKaonTruthTree::AddNeutralKaonDaughter2TruthVariables(output(), 10);
  // Function in standardPDTree.cxx where the truth tree variables are defined: momentum, pdg, etc.
  // Function in standardPDTree.cxx -> beamParticleTruthDaughters()
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
  // Fill truth tree for all particles to include vertex information
  if (!part) return false;

  if(part->PDG == 310){
    // std::cout << "DEBUG: Filling truth tree for K0" << std::endl;
    return true;
  }
  else{
    return false;
  }
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************
    // Fill the common variables
    pdBaseAnalysis::FillTruthTree(part);

    // If this is a true K0S/K0L/K0 decaying into 2 pi0, print run/subrun/event
    const int absPdg = std::abs(part.PDG);
    const bool isK0S = (absPdg == 310);
    const bool isK0L = (absPdg == 130);
    const bool isK0  = (absPdg == 311);

    if ((isK0S || isK0L || isK0) &&
        part.ProcessEnd == AnaTrueParticleB::Decay &&
        part.Daughters.size() == 2) {

      AnaTrueParticlePD* dau1 = nullptr;
      AnaTrueParticlePD* dau2 = nullptr;

      for (int i = 0; i < GetSpill().TrueParticles.size(); ++i) {
        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[i]);
        if (!truePart) continue;
        if (truePart->ID == part.Daughters[0]) dau1 = truePart;
        if (truePart->ID == part.Daughters[1]) dau2 = truePart;
      }

      if (dau1 && dau2 && dau1->PDG == 111 && dau2->PDG == 111) {
        const char* k0Type = isK0S ? "K0S" : (isK0L ? "K0L" : "K0");
        const AnaEventInfoPD* evtInfo = static_cast<const AnaEventInfoPD*>(GetEvent().EventInfo);
        if (evtInfo) {
          std::cout << "[neutralKaonAnalysis] True " << k0Type << " -> pi0 pi0 decay in event: "
                    << "Run " << evtInfo->Run
                    << " SubRun " << evtInfo->SubRun
                    << " Event " << evtInfo->Event
                    << std::endl;
        } else {
          std::cout << "[neutralKaonAnalysis] True " << k0Type << " -> pi0 pi0 decay in current event" << std::endl;
        }
      }
    }

    // The truth tree is meant for individual particle information, not analysis results
    // Vertex candidates are analysis results and belong in the ana tree only

    // OPTIMIZATION: Build hash maps once for O(1) lookups
    std::unordered_map<Int_t, AnaTrueParticlePD*> trueParticleByID;
    std::unordered_map<Int_t, AnaParticlePD*> recoParticleByTrueID;

    // Build true particle map
    for(int i = 0; i < GetSpill().TrueParticles.size(); i++){
      AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(GetSpill().TrueParticles[i]);
      if(truePart){
        trueParticleByID[truePart->ID] = truePart;
      }
    }

    // Build reco particle map (by associated true ID)
    for(size_t i = 0; i < GetBunch().Particles.size(); i++){
      AnaParticlePD* recoPart = static_cast<AnaParticlePD*>(GetBunch().Particles[i]);
      if(recoPart && recoPart->TrueObject){
        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(recoPart->TrueObject);
        if(truePart){
          recoParticleByTrueID[truePart->ID] = recoPart;
        }
      }
    }

    // Check if K0 has a reconstructed object - O(1) hash map lookup
    AnaParticlePD* k0RecoObject = nullptr;
    auto it_k0 = recoParticleByTrueID.find(part.ID);
    if(it_k0 != recoParticleByTrueID.end()){
      k0RecoObject = it_k0->second;
    }
    neutralKaonTruthTree::FillNeutralKaonTruthVariables(output(), part, k0RecoObject != nullptr);

    // Fill vertex reconstruction debugging variables
    // Get reco particles for K0 daughters (if they exist) - O(1) hash map lookups
    AnaParticlePD* daughter1Reco = nullptr;
    AnaParticlePD* daughter2Reco = nullptr;
    AnaParticlePD* parentReco = nullptr;

    if(part.Daughters.size() > 0){
      // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
      auto it_d1 = trueParticleByID.find(part.Daughters[0]);
      if(it_d1 != trueParticleByID.end()){
        AnaTrueParticlePD* daughter1True = it_d1->second;
        if(daughter1True){
          auto it_d1r = recoParticleByTrueID.find(daughter1True->ID);
          if(it_d1r != recoParticleByTrueID.end()){
            daughter1Reco = it_d1r->second;
          }
        }
      }

      if(part.Daughters.size() > 1){
        // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
        auto it_d2 = trueParticleByID.find(part.Daughters[1]);
        if(it_d2 != trueParticleByID.end()){
          AnaTrueParticlePD* daughter2True = it_d2->second;
          if(daughter2True){
            auto it_d2r = recoParticleByTrueID.find(daughter2True->ID);
            if(it_d2r != recoParticleByTrueID.end()){
              daughter2Reco = it_d2r->second;
            }
          }
        }
      }
    }

    // Get parent reco particle if it exists - O(1) hash map lookup
    auto it_par = trueParticleByID.find(part.ParentID);
    if(it_par != trueParticleByID.end()){
      AnaTrueParticlePD* parentTrue = it_par->second;
      if(parentTrue){
        auto it_parr = recoParticleByTrueID.find(parentTrue->ID);
        if(it_parr != recoParticleByTrueID.end()){
          parentReco = it_parr->second;
        }

        // If not found in regular particles, check if it's the beam particle
        if(!parentReco){
          AnaParticlePD* beamParticle = static_cast<AnaParticlePD*>(static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle);
          if(beamParticle && beamParticle->TrueObject){
            AnaTrueParticlePD* trueBeam = static_cast<AnaTrueParticlePD*>(beamParticle->TrueObject);
            if(trueBeam && trueBeam->ID == parentTrue->ID){
              parentReco = beamParticle;
            }
          }
        }
      }
    }

    // Get parameters for vertex reconstruction
    double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
    double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

    // Fill the debugging variables
    neutralKaonTruthTree::FillVertexReconstructionDebugVariables(output(), part, daughter1Reco, daughter2Reco,
                                                                 parentReco, maxDaughterDistance, trackFitLength, maxDaughterDistance);

    // Fill parent information if it exists - O(1) hash map lookup
    auto it_parent = trueParticleByID.find(part.ParentID);
    if(it_parent != trueParticleByID.end()){
      AnaTrueParticlePD* parent = it_parent->second;
      if(parent){
        // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
        AnaParticlePD* parentRecoObject = nullptr;
        auto it_parentr = recoParticleByTrueID.find(parent->ID);
        if(it_parentr != recoParticleByTrueID.end()){
          parentRecoObject = it_parentr->second;
        }

        // If not found in regular particles, check if it's the beam particle
        if(!parentRecoObject){
          AnaParticlePD* beamParticle = static_cast<AnaParticlePD*>(static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle);
          if(beamParticle && beamParticle->TrueObject){
            AnaTrueParticlePD* trueBeam = static_cast<AnaTrueParticlePD*>(beamParticle->TrueObject);
            if(trueBeam && trueBeam->ID == parent->ID){
              parentRecoObject = beamParticle;
            }
          }
        }
        neutralKaonTruthTree::FillNeutralKaonParentTruthVariables(output(), *parent, parentRecoObject != nullptr);
      }
      else{
        // Parent not found - fill hasrecoobject with 0
        output().FillVectorVar(neutralKaonTruthTree::k0parhasrecoobject, 0);
      }
    }
    else{
      // Parent not found - fill hasrecoobject with 0
      output().FillVectorVar(neutralKaonTruthTree::k0parhasrecoobject, 0);
    }

    // Fill daughter information if daughters exist - O(1) hash map lookups
    if(part.Daughters.size() > 0){
      // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
      auto it_dau1 = trueParticleByID.find(part.Daughters[0]);
      if(it_dau1 != trueParticleByID.end()){
        AnaTrueParticlePD* daughter1 = it_dau1->second;
        if(daughter1){
          AnaParticlePD* daughter1RecoObject = nullptr;
          auto it_dau1r = recoParticleByTrueID.find(daughter1->ID);
          if(it_dau1r != recoParticleByTrueID.end()){
            daughter1RecoObject = it_dau1r->second;
          }
          neutralKaonTruthTree::FillNeutralKaonDaughter1TruthVariables(output(), *daughter1, daughter1RecoObject != nullptr);
        }
        else{
          // Daughter1 not found - fill hasrecoobject with 0
          output().FillVectorVar(neutralKaonTruthTree::k0dau1hasrecoobject, 0);
        }
      }
      else{
        // Daughter1 not found - fill hasrecoobject with 0
        output().FillVectorVar(neutralKaonTruthTree::k0dau1hasrecoobject, 0);
      }

      if(part.Daughters.size() > 1){
        // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
        auto it_dau2 = trueParticleByID.find(part.Daughters[1]);
        if(it_dau2 != trueParticleByID.end()){
          AnaTrueParticlePD* daughter2 = it_dau2->second;
          if(daughter2){
            AnaParticlePD* daughter2RecoObject = nullptr;
            auto it_dau2r = recoParticleByTrueID.find(daughter2->ID);
            if(it_dau2r != recoParticleByTrueID.end()){
              daughter2RecoObject = it_dau2r->second;
            }
            neutralKaonTruthTree::FillNeutralKaonDaughter2TruthVariables(output(), *daughter2, daughter2RecoObject != nullptr);
          }
          else{
            // Daughter2 not found - fill hasrecoobject with 0
            output().FillVectorVar(neutralKaonTruthTree::k0dau2hasrecoobject, 0);
          }
        }
        else{
          // Daughter2 not found - fill hasrecoobject with 0
          output().FillVectorVar(neutralKaonTruthTree::k0dau2hasrecoobject, 0);
        }
      }
      else{
        // No daughter2 - fill hasrecoobject with 0
        output().FillVectorVar(neutralKaonTruthTree::k0dau2hasrecoobject, 0);
      }
    }
    else{
      // No daughters at all - fill both hasrecoobject with 0
      output().FillVectorVar(neutralKaonTruthTree::k0dau1hasrecoobject, 0);
      output().FillVectorVar(neutralKaonTruthTree::k0dau2hasrecoobject, 0);
    }

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
    }
  }

}