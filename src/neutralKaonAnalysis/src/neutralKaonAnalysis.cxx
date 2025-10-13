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
  // standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  // standardPDTree::AddStandardVariables_BeamParticleReco(output());
  // standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
  // standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),50);
  // standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),50);
  // standardPDTree::AddStandardVariables_BeamTruthDaughters(output(),50);

  // AddVarI(output(), seltrk_dau_trueparentpdg, "Parent PDG of reco daughter");

  // AddVarI(output(), nAllParticles, "Number of all particles with valid start positions");

  // Add neutral particle candidates variables
  // neutralKaonTree::AddNeutralKaonVariables_Candidates(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0Par(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0Vtx(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0Brother(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0vtxDaughter1(output(), 1000);
  neutralKaonTree::AddNeutralKaonVariables_K0vtxDaughter2(output(), 1000);

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

  // ---------- Additional candidate variables --------------
  // if(box().MainTrack){
  //   // Fill beam particle information
  //   standardPDTree::FillStandardVariables_BeamParticleReco(output(), box().MainTrack);
  //   standardPDTree::FillStandardVariables_BeamParticleTrue(output(), box().MainTrack);
  //   standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);

  //   // Fill beam particle daughters information
  //   int ndau = std::min(50, (int)box().MainTrack->Daughters.size());
  //   for(int i = 0; i < ndau; i++){
  //     AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]);
  //     if (!daughter) continue;

  //     standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), daughter);
  //     standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), daughter);

  //     // Fill parent PDG information
  //     if (daughter->TrueObject) {
  //       AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(daughter->TrueObject);
  //       if (trueDaughter) {
  //         output().FillVar(seltrk_dau_trueparentpdg, trueDaughter->ParentPDG);
  //       }
  //     }

  //     output().IncrementCounter(standardPDTree::seltrk_ndau);
  //   }

  //   // Fill truth daughter counter (needed by other variables)
  //   AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
  //   if (trueBeamPart) {
  //     int ndau_truth = (int)trueBeamPart->Daughters.size();
  //     for(int j = 0; j < ndau_truth; j++){
  //       AnaTrueParticlePD* truthdau = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, trueBeamPart->Daughters[j]);
  //       if (!truthdau) continue;
  //       output().IncrementCounter(standardPDTree::seltrk_truthdau_ndau);
  //     }
  //   }
  // }

  // Fill vertex candidates information for all events (regardless of main track)
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());
  // output().FillVar(nAllParticles, neutralKaonBox.nAllParticles);

  // Fill neutral particle candidates data
  // neutralKaonTree::FillNeutralKaonVariables_Candidates(output(), neutralKaonBox.neutralParticleCandidates, GetEvent());

  // Fill the number of neutral particle candidates
  // output().FillVar(neutralKaonTree::nk0, static_cast<Int_t>(neutralKaonBox.neutralParticleCandidates.size()));

  // Fill individual candidate data
  if(neutralKaonBox.neutralParticleCandidates.size() > 0){
    for(size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++){
      neutralKaonTree::FillNeutralKaonVariables(output(), neutralKaonBox.neutralParticleCandidates[i], GetEvent(), GetSpill().Beam);
      output().IncrementCounter(neutralKaonTree::nk0);
    }
  }

  // Create event display for events that pass the selection and have neutral particles with true objects
  if (_EventDisplay) {
    // Check if any neutral particle candidate has an associated true particle
    bool hasNeutralWithTrue = false;
    const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());

    for (size_t i = 0; i < neutralKaonBox.neutralParticleCandidates.size(); i++) {
      // Check if this neutral particle has an associated true particle that is K0, Pi0, or Gamma
      AnaNeutralParticlePD* neutralParticle = neutralKaonBox.neutralParticleCandidates[i];
      if (neutralParticle->TrueObject != nullptr) {
        AnaTrueParticleB* trueNeutralParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        // Check for K0 (310), Pi0 (111), or Gamma (22)
        if (trueNeutralParticle->PDG == 310 || trueNeutralParticle->PDG == 130) {
          if (trueNeutralParticle->ProcessEnd == 2) {
            hasNeutralWithTrue = true;
            break;
          }
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

    // The truth tree is meant for individual particle information, not analysis results
    // Vertex candidates are analysis results and belong in the ana tree only

    // Check if K0 has a reconstructed object
    AnaParticlePD* k0RecoObject = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, part.ID);
    neutralKaonTruthTree::FillNeutralKaonTruthVariables(output(), part, k0RecoObject != nullptr);

    // Fill vertex reconstruction debugging variables
    // Get reco particles for K0 daughters (if they exist)
    AnaParticlePD* daughter1Reco = nullptr;
    AnaParticlePD* daughter2Reco = nullptr;
    AnaParticlePD* parentReco = nullptr;

    if(part.Daughters.size() > 0){
      AnaTrueParticlePD* daughter1True = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.Daughters[0]);
      if(daughter1True){
        daughter1Reco = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, daughter1True->ID);
      }

      if(part.Daughters.size() > 1){
        AnaTrueParticlePD* daughter2True = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.Daughters[1]);
        if(daughter2True){
          daughter2Reco = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, daughter2True->ID);
        }
      }
    }

    // Get parent reco particle if it exists
    AnaTrueParticlePD* parentTrue = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.ParentID);
    if(parentTrue){
      parentReco = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, parentTrue->ID);

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

    // Get parameters for vertex reconstruction
    double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
    double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
    double vertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");

    // Fill the debugging variables
    neutralKaonTruthTree::FillVertexReconstructionDebugVariables(output(), part, daughter1Reco, daughter2Reco,
                                                                 parentReco, maxDaughterDistance, trackFitLength, vertexRadius);

    // Fill parent information if it exists
    AnaTrueParticlePD* parent = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.ParentID);
    if(parent){
      AnaParticlePD* parentRecoObject = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, parent->ID);

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

    // Fill daughter information if daughters exist
    if(part.Daughters.size() > 0){
      AnaTrueParticlePD* daughter1 = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.Daughters[0]);
      if(daughter1){
        AnaParticlePD* daughter1RecoObject = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, daughter1->ID);
        neutralKaonTruthTree::FillNeutralKaonDaughter1TruthVariables(output(), *daughter1, daughter1RecoObject != nullptr);
      }
      else{
        // Daughter1 not found - fill hasrecoobject with 0
        output().FillVectorVar(neutralKaonTruthTree::k0dau1hasrecoobject, 0);
      }

      if(part.Daughters.size() > 1){
        AnaTrueParticlePD* daughter2 = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, part.Daughters[1]);
        if(daughter2){
          AnaParticlePD* daughter2RecoObject = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(GetBunch().Particles, daughter2->ID);
          neutralKaonTruthTree::FillNeutralKaonDaughter2TruthVariables(output(), *daughter2, daughter2RecoObject != nullptr);
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
      neutralKaonAnaUtils::FillNeutralParticleSignalBackgroundCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonAnaUtils::FillSignalCandidateCategory(neutralKaonBox.neutralParticleCandidates[i]);
      neutralKaonAnaUtils::FillNeutralParticlePDGCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonAnaUtils::FillNeutralParticleChargeCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());
      neutralKaonAnaUtils::FillNeutralParticleNoTruthCategory(neutralKaonBox.neutralParticleCandidates[i], GetEvent());

    }
  }

}