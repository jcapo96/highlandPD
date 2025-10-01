#include "neutralKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "neutralKaonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"
#include "ToyBoxNeutralKaon.hxx"

#include "pdAnalysisUtils.hxx"
#include "standardPDTree.hxx"
#include <sstream>
#include <iostream>

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"

#include "baseToyMaker.hxx"

//********************************************************************
neutralKaonAnalysis::neutralKaonAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************

  _EventDisplay = nullptr;

  // Add the package version
  //  ND::versioning().AddPackage("StoppingProtonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
}

//********************************************************************
neutralKaonAnalysis::~neutralKaonAnalysis() {
//********************************************************************
  if (_EventDisplay) {
    delete _EventDisplay;
    _EventDisplay = nullptr;
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

  // Initialize event display
  _EventDisplay = new pdEventDisplay();

  // Get event display parameters
  bool createEventDisplay = ND::params().GetParameterI("neutralKaonAnalysis.CreateEventDisplay");
  bool saveToRootFile = ND::params().GetParameterI("neutralKaonAnalysis.SaveToRootFile");
  std::string outputDirectory = ND::params().GetParameterS("neutralKaonAnalysis.OutputDirectory");
  int maxEventsToDisplay = ND::params().GetParameterI("neutralKaonAnalysis.MaxEventsToDisplay");
  double eventPercentage = ND::params().GetParameterD("neutralKaonAnalysis.EventDisplayPercentage");

  // Parse required particle PDGs
  std::vector<int> requiredPDGs;
  std::string pdgString = ND::params().GetParameterS("neutralKaonAnalysis.RequiredParticlePDGs");
  if (!pdgString.empty()) {
    std::istringstream iss(pdgString);
    std::string token;
    while (std::getline(iss, token, ',')) {
      try {
        int pdg = std::stoi(token);
        requiredPDGs.push_back(pdg);
      } catch (const std::exception& e) {
        std::cout << "Warning: Could not parse PDG code: " << token << std::endl;
      }
    }
  }

  double vertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");
  int minVertexDaughters = ND::params().GetParameterI("neutralKaonAnalysis.MinVertexDaughters");

  // Initialize the event display
  _EventDisplay->Initialize(createEventDisplay, saveToRootFile, outputDirectory,
                           maxEventsToDisplay, eventPercentage, requiredPDGs,
                           vertexRadius, minVertexDaughters);

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("bestcandidate");

  // Add custom categories for neutral particle analysis
  neutralKaonAnaUtils::AddCustomCategories();

  return true;
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
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());
  standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
  standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),50);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),50);
  standardPDTree::AddStandardVariables_BeamTruthDaughters(output(),50);

  AddVarI(output(), seltrk_dau_trueparentpdg, "Parent PDG of reco daughter");

  AddVarI(output(), nAllParticles, "Number of all particles with valid start positions");

  // Add neutral particle candidates variables
  neutralKaonTree::AddNeutralKaonVariables_Candidates(output(), 1000);

}

//********************************************************************
void neutralKaonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  pdBaseAnalysis::DefineTruthTree();
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

  // ---------- Additional candidate variables --------------
  if(box().MainTrack){
    // Fill beam particle information
    standardPDTree::FillStandardVariables_BeamParticleReco(output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleTrue(output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);

    // Fill beam particle daughters information
    int ndau = std::min(50, (int)box().MainTrack->Daughters.size());
    for(int i = 0; i < ndau; i++){
      AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]);
      if (!daughter) continue;

      standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), daughter);
      standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), daughter);

      // Fill parent PDG information
      if (daughter->TrueObject) {
        AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(daughter->TrueObject);
        if (trueDaughter) {
          output().FillVar(seltrk_dau_trueparentpdg, trueDaughter->ParentPDG);
        }
      }

      output().IncrementCounter(standardPDTree::seltrk_ndau);
    }

    // Fill truth daughter counter (needed by other variables)
    AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
    if (trueBeamPart) {
      int ndau_truth = (int)trueBeamPart->Daughters.size();
      for(int j = 0; j < ndau_truth; j++){
        AnaTrueParticlePD* truthdau = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, trueBeamPart->Daughters[j]);
        if (!truthdau) continue;
        output().IncrementCounter(standardPDTree::seltrk_truthdau_ndau);
      }
    }
  }

  // Fill vertex candidates information for all events (regardless of main track)
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());
  output().FillVar(nAllParticles, neutralKaonBox.nAllParticles);

  // Fill neutral particle candidates data
  neutralKaonTree::FillNeutralKaonVariables_Candidates(output(), neutralKaonBox.neutralParticleCandidates, GetEvent());

  // Fill individual candidate data
  if(neutralKaonBox.neutralParticleCandidates.size() > 0){
    for(size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++){
      neutralKaonTree::FillNeutralKaonVariables_SingleCandidate(output(), neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonTree::FillNeutralKaonVariables_Candidates(output(), neutralKaonBox.neutralParticleCandidates, GetEvent());
      output().IncrementCounter(neutralKaonTree::nk0);
    }
  }

  // Create event display for events that pass the selection and have neutral particles with true objects
  if (_EventDisplay && _EventDisplay->ShouldCreateEventDisplay()) {
    // Check if any neutral particle candidate has an associated true particle
    bool hasNeutralWithTrue = false;
    const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());

    for (size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++) {
      // Check if this neutral particle has an associated true particle that is K0, Pi0, or Gamma
      AnaNeutralParticlePD* neutralParticle = neutralKaonBox.neutralParticleCandidates[i];
      if (neutralParticle->TrueObject != nullptr) {
        AnaTrueParticleB* trueNeutralParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        // Check for K0 (310), Pi0 (111), or Gamma (22)
        if (trueNeutralParticle->PDG == 310 || trueNeutralParticle->PDG == 111 || trueNeutralParticle->PDG == 22) {
          hasNeutralWithTrue = true;
          break;
        }
      }
    }

    // Only create event display if there's at least one neutral particle with a true object
    if (hasNeutralWithTrue) {
      int eventNumber = GetEvent().EventInfo->Event;
      _EventDisplay->CreateEventDisplay(GetEvent(), eventNumber, neutralKaonBox);
    }
  }

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

  return true; // Fill truth tree for all particles
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************
    // Fill the common variables
    pdBaseAnalysis::FillTruthTree(part);

    // The truth tree is meant for individual particle information, not analysis results
    // Vertex candidates are analysis results and belong in the ana tree only
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
      neutralKaonAnaUtils::FillNeutralParticleSignalBackgroundCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonAnaUtils::FillNeutralParticleDetailedBackgroundCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
    }
  }

}