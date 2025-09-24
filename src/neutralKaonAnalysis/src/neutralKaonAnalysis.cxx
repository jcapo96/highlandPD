#include "neutralKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "neutralKaonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"

#include "pdAnalysisUtils.hxx"
#include "standardPDTree.hxx"
#include <sstream>
#include <iostream>

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

#include "ParticlePositionSCECorrection.hxx"
#include "SCEGeometricVariation.hxx"
#include "pdEventDisplay.hxx"

#include "baseToyMaker.hxx"

//********************************************************************
neutralKaonAnalysis::neutralKaonAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************

  // Initialize event display pointer
  _eventDisplay = NULL;

  // Add the package version
  //  ND::versioning().AddPackage("StoppingProtonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
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

    // Event display parameters
    _CreateEventDisplay = ND::params().GetParameterI("neutralKaonAnalysis.CreateEventDisplay");
    _SaveToRootFile = ND::params().GetParameterI("neutralKaonAnalysis.SaveToRootFile");
    _OutputDirectory = ND::params().GetParameterS("neutralKaonAnalysis.OutputDirectory");
    _MaxEventsToDisplay = ND::params().GetParameterI("neutralKaonAnalysis.MaxEventsToDisplay");
    _EventDisplayPercentage = ND::params().GetParameterD("neutralKaonAnalysis.EventDisplayPercentage");
    _VertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.VertexRadius");
    _MinVertexDaughters = ND::params().GetParameterI("neutralKaonAnalysis.MinVertexDaughters");
    _OnlySignalEvents = ND::params().GetParameterI("neutralKaonAnalysis.OnlySignalEvents");

    // Parse required particle PDGs from parameters
    std::string requiredPDGsStr = ND::params().GetParameterS("neutralKaonAnalysis.RequiredParticlePDGs");
    _RequiredParticlePDGs.clear();
    if (!requiredPDGsStr.empty() && requiredPDGsStr != "none") {
        // Remove quotes if present
        if (requiredPDGsStr.front() == '"' && requiredPDGsStr.back() == '"') {
            requiredPDGsStr = requiredPDGsStr.substr(1, requiredPDGsStr.length() - 2);
        }

        std::istringstream iss(requiredPDGsStr);
        std::string pdgStr;
        while (std::getline(iss, pdgStr, ',')) {
            // Trim whitespace from the string
            pdgStr.erase(0, pdgStr.find_first_not_of(" \t"));
            pdgStr.erase(pdgStr.find_last_not_of(" \t") + 1);

            // Only parse if the string is not empty
            if (!pdgStr.empty()) {
                try {
                    int pdg = std::stoi(pdgStr);
                    _RequiredParticlePDGs.push_back(pdg);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Warning: Could not parse PDG '" << pdgStr << "' from RequiredParticlePDGs parameter. Skipping." << std::endl;
                }
            }
        }
    }

    // Debug: Print loaded PDG requirements
    if (!_RequiredParticlePDGs.empty()) {
        std::cout << "Event display will only save events containing these particle types: ";
        for (size_t i = 0; i < _RequiredParticlePDGs.size(); ++i) {
            std::cout << _RequiredParticlePDGs[i];
            if (i < _RequiredParticlePDGs.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Event display will save events with any particle types (no filtering)" << std::endl;
    }

    // Initialize event display
    _eventDisplay = new pdEventDisplay();
    _eventDisplay->Initialize(_CreateEventDisplay, _SaveToRootFile, _OutputDirectory, _MaxEventsToDisplay, _EventDisplayPercentage, _RequiredParticlePDGs, _VertexRadius, _MinVertexDaughters);

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("bestcandidate");

  // Add custom categories for neutral kaon analysis
  AddVertexPionPairCategory();
  AddVertexParticleCountCategory();
  AddK0InVtxCategory();
  AddVertexParentPDGCategory();

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

  // // Add vertex candidates variables
  // // Increased from 500 to 1000 to accommodate the new combinatorial vertex creation algorithm
  neutralKaonTree::AddNeutralKaonVariables_VertexCandidates(output(), 1000);

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

  // // Initialize counters for this event
  // output().InitializeCounter(neutralKaonTree::n_recovtx_candidates);
  // output().InitializeCounter(neutralKaonTree::n_true_vertex_candidates);

  // Fill vertex candidates data
  // const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());
  neutralKaonTree::FillNeutralKaonVariables_VertexCandidates(output(), neutralKaonBox.reconVertexCandidates, neutralKaonBox.trueVertexCandidates);

  // Create event display if enabled and event passes the minimum accumulated level to save
  if (_eventDisplay && _eventDisplay->ShouldCreateEventDisplay()) {
    // Only create event display for events that pass the minimum accumulated level to save
    SelectionBase* neutralKaonSelection = sel().GetSelection("neutralKaonSelection");
    if (neutralKaonSelection) {
      int maxAccumLevel = 0;
      for(UInt_t ibranch = 0; ibranch < neutralKaonSelection->GetNBranches(); ibranch++){
        if(neutralKaonSelection->GetAccumCutLevel(ibranch) > maxAccumLevel){
          maxAccumLevel = neutralKaonSelection->GetAccumCutLevel(ibranch);
        }
      }

      // Get the minimum accumulated level to save from parameters
      int minAccumLevelToSave = ND::params().GetParameterI("neutralKaonAnalysis.MinAccumLevelToSave");

      // Only display events that pass the minimum accumulated level to save
      if (maxAccumLevel >= minAccumLevelToSave) {
        // Check if event contains required particles BEFORE creating event display
        _eventDisplay->CreateEventDisplay(GetEvent(), GetEvent().EventInfo->Event);
      }
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

  // Fill vertex pion pair categories for reconstructed vertex candidates
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());

  // Fill categories for the best reconstructed vertex candidate if available
  // Fill categories for all vertices
  FillK0InVtxCategory(neutralKaonBox.reconVertexCandidates);

  // Fill other categories for best vertex if available
  // if(neutralKaonBox.BestReconVertexCandidateIndex >= 0 &&
  //    neutralKaonBox.BestReconVertexCandidateIndex < (int)neutralKaonBox.reconVertexCandidates.size()) {
  //   AnaVertexPD* bestVertex = neutralKaonBox.reconVertexCandidates[neutralKaonBox.BestReconVertexCandidateIndex];
  //   if(bestVertex) {
  //     FillVertexPionPairCategory(bestVertex);
  //     FillVertexParticleCountCategory(bestVertex);
  //     FillVertexParentPDGCategory(bestVertex);
  //   }
  // }
  // If no best vertex, try to fill for any available vertex
  // else if(!neutralKaonBox.reconVertexCandidates.empty()) {
  //   AnaVertexPD* firstVertex = neutralKaonBox.reconVertexCandidates[0];
  //   if(firstVertex) {
  //     FillVertexPionPairCategory(firstVertex);
  //     FillVertexParticleCountCategory(firstVertex);
  //     FillVertexParentPDGCategory(firstVertex);
  //   }
  // }
  // If no vertices available, set to no truth
  // else {
  //   FillVertexPionPairCategory(nullptr);
  //   FillVertexParticleCountCategory(nullptr);
  //   FillVertexParentPDGCategory(nullptr);
  // }
}

//********************************************************************
void neutralKaonAnalysis::AddVertexPionPairCategory(){
//********************************************************************

  std::string part_types[] = {"signal", "background", NAMEOTHER};
  int part_codes[]         = {1        , 0           , CATOTHER};
  int part_colors[]        = {2        , 4           , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddCategory("vertexpionpair", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void neutralKaonAnalysis::AddVertexParticleCountCategory(){
//********************************************************************

  std::string part_types[] = {"0_daughters", "1_daughter", "2_daughters", "3_daughters", "4_daughters", "more_than_4", NAMEOTHER};
  int part_codes[]         = {0            , 1           , 2            , 3            , 4            , 5            , CATOTHER};
  int part_colors[]        = {1            , 2           , 3            , 4            , 6            , 7            , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddCategory("vertexparticlecount", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void neutralKaonAnalysis::AddK0InVtxCategory(){
//********************************************************************

  std::string part_types[] = {"signal", "background", NAMEOTHER};
  int part_codes[]         = {1        , 0           , CATOTHER};
  int part_colors[]        = {3        , 2           , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("k0invtx", neutralKaonTree::n_recovtx_candidates, "n_recovtx_candidates",
    NPART, part_types, part_codes, part_colors,
    1, -1000);
  // anaUtils::_categ->AddObjectCategory("k0invtx", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void neutralKaonAnalysis::FillVertexPionPairCategory(AnaVertexPD* vertex){
//********************************************************************

  if(!vertex) {
    anaUtils::_categ->SetCode("vertexpionpair", CATNOTRUTH, CATNOTRUTH);
    return;
  }

  // Check if there are exactly two pions (pi+ and pi-) in the reconstructed particles
  bool hasPiPlus = false;
  bool hasPiMinus = false;

  // Loop through all reconstructed particles in the vertex
  for(int i = 0; i < vertex->NParticles; i++) {
    AnaParticlePD* recoPart = vertex->Particles[i];
    if(!recoPart) continue;

    // Get the associated true particle
    AnaTrueParticleB* truePart = recoPart->GetTrueParticle();
    if(!truePart) continue;

    if(truePart->PDG == 211) {      // pi+
      hasPiPlus = true;
    }
    else if(truePart->PDG == -211) { // pi-
      hasPiMinus = true;
    }
  }

  // Set category based on pion pair presence
  if(hasPiPlus && hasPiMinus) {
    anaUtils::_categ->SetCode("vertexpionpair", 1, CATOTHER); // signal
  }
  else {
    anaUtils::_categ->SetCode("vertexpionpair", 0, CATOTHER); // background
  }
}

//********************************************************************
void neutralKaonAnalysis::FillVertexParticleCountCategory(AnaVertexPD* vertex){
//********************************************************************

  if(!vertex) {
    anaUtils::_categ->SetCode("vertexparticlecount", CATNOTRUTH, CATNOTRUTH);
    return;
  }

  // Count the number of daughters in the vertex (excluding the parent)
  // The first particle is the parent, so daughters = total particles - 1
  int daughterCount = vertex->NParticles - 1;

  // Set category based on daughter count
  if(daughterCount == 0) {
    anaUtils::_categ->SetCode("vertexparticlecount", 0, CATOTHER); // 0 daughters
  }
  else if(daughterCount == 1) {
    anaUtils::_categ->SetCode("vertexparticlecount", 1, CATOTHER); // 1 daughter
  }
  else if(daughterCount == 2) {
    anaUtils::_categ->SetCode("vertexparticlecount", 2, CATOTHER); // 2 daughters
  }
  else if(daughterCount == 3) {
    anaUtils::_categ->SetCode("vertexparticlecount", 3, CATOTHER); // 3 daughters
  }
  else if(daughterCount == 4) {
    anaUtils::_categ->SetCode("vertexparticlecount", 4, CATOTHER); // 4 daughters
  }
  else if(daughterCount > 4) {
    anaUtils::_categ->SetCode("vertexparticlecount", 5, CATOTHER); // more than 4 daughters
  }
  else {
    // This should not happen (negative daughter count), but handle it
    anaUtils::_categ->SetCode("vertexparticlecount", CATOTHER, CATOTHER);
  }
}

//********************************************************************
void neutralKaonAnalysis::FillK0InVtxCategory(const std::vector<AnaVertexPD*>& vertices){
//********************************************************************

  // Loop over all vertices and classify each one
  for(const auto& vertex : vertices) {
    if(!vertex) {
      anaUtils::_categ->SetObjectCode("k0invtx", CATNOTRUTH, CATNOTRUTH, -1);
      continue;
    }

    // Check if the two daughters are π+ and π- from the same K0 (PDG=310) by checking both PDG and particle ID
    bool isSignal = false;

    // We need exactly 2 daughters (plus 1 parent = 3 total particles)
    if(vertex->NParticles >= 3) {
      std::vector<AnaTrueParticleB*> daughterTrueParts;

      // Skip the first particle (parent) and collect the daughters
      for(int i = 1; i < vertex->NParticles; i++) {
        AnaParticlePD* recoPart = vertex->Particles[i];
        if(!recoPart) continue;

        // Get the associated true particle
        AnaTrueParticleB* truePart = recoPart->GetTrueParticle();
        if(!truePart) continue;

        daughterTrueParts.push_back(truePart);
      }

      // Check if we have exactly 2 daughters
      if(daughterTrueParts.size() == 2) {
        AnaTrueParticleB* daughter1 = daughterTrueParts[0];
        AnaTrueParticleB* daughter2 = daughterTrueParts[1];

        // Check if both daughters have the same K0 parent (PDG=310 and same ParentID)
        if(daughter1->ParentPDG == 310 && daughter2->ParentPDG == 310) {
          // Check if they have the same parent particle ID (same K0)
          if(daughter1->ParentID == daughter2->ParentID && daughter1->ParentID != 0) {
            // Now check if they are specifically π+ and π- (signal definition)
            if((daughter1->PDG == 211 && daughter2->PDG == -211) ||
               (daughter1->PDG == -211 && daughter2->PDG == 211)) {

              // Additional condition: Check if the K0 parent has a K+ parent (PDG=321)
              // and that this K+ is the same particle that originated the vertex

              // Find the K0 parent particle in the true particles list
              AnaTrueParticleB* k0Parent = nullptr;
              AnaTrueParticleB** allTrueParticles = GetEvent().TrueParticles;
              int nTrueParticles = GetEvent().nTrueParticles;

              for (int j = 0; j < nTrueParticles; j++) {
                if (allTrueParticles[j] && allTrueParticles[j]->ID == daughter1->ParentID) {
                  k0Parent = allTrueParticles[j];
                  break;
                }
              }

              if (k0Parent) {
                // Check if the K0 parent has a K+ parent (PDG=321)
                if (k0Parent->ParentPDG == 321) {
                  // Find the K+ parent particle
                  AnaTrueParticleB* kPlusParent = nullptr;
                  for (int k = 0; k < nTrueParticles; k++) {
                    if (allTrueParticles[k] && allTrueParticles[k]->ID == k0Parent->ParentID) {
                      kPlusParent = allTrueParticles[k];
                      break;
                    }
                  }

                  if (kPlusParent) {
                    // Check if this K+ is the same particle that originated the vertex
                    // The vertex parent should be the reconstructed particle corresponding to this K+
                    AnaParticlePD* vertexParent = static_cast<AnaParticlePD*>(vertex->Particles[0]);
                    if (vertexParent) {
                      AnaTrueParticleB* vertexParentTrue = vertexParent->GetTrueParticle();
                      if (vertexParentTrue && vertexParentTrue->ID == kPlusParent->ID) {
                        isSignal = true; // All conditions met for signal vertex
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Set category based on signal/background classification for this vertex
    if(isSignal) {
      anaUtils::_categ->SetObjectCode("k0invtx", 1, CATOTHER, -1); // signal: daughters from same K0
    }
    else {
      anaUtils::_categ->SetObjectCode("k0invtx", 0, CATOTHER, -1); // background: daughters not from same K0
    }
  }
}

//********************************************************************
void neutralKaonAnalysis::AddVertexParentPDGCategory(){
//********************************************************************

  std::string part_types[] = {"K+", "K-", "K0", "pi+", "pi-", "pi0", "p", "n", "mu+", "mu-", "e+", "e-", "gamma", "other", NAMEOTHER};
  int part_codes[]         = {321  , -321 , 310 , 211  , -211  , 111  , 2212, 2112, 13   , -13   , 11  , -11  , 22    , 0     , CATOTHER};
  int part_colors[]        = {2    , 3   , 4   , 5    , 6    , 7    , 8   , 9   , 10   , 11   , 12  , 13  , 14    , 15    , COLOTHER};

  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  anaUtils::_categ->AddCategory("vertexparentpdg", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void neutralKaonAnalysis::FillVertexParentPDGCategory(AnaVertexPD* vertex){
//********************************************************************

  if(!vertex) {
    anaUtils::_categ->SetCode("vertexparentpdg", CATNOTRUTH, CATNOTRUTH);
    return;
  }

  // Get the parent particle of the vertex
  AnaParticlePD* parent = vertex->Parent;
  if(!parent) {
    anaUtils::_categ->SetCode("vertexparentpdg", CATNOTRUTH, CATNOTRUTH);
    return;
  }

  // Get the true particle associated with the parent
  AnaTrueParticleB* trueParent = parent->GetTrueParticle();
  if(!trueParent) {
    anaUtils::_categ->SetCode("vertexparentpdg", CATNOTRUTH, CATNOTRUTH);
    return;
  }

  // Set the category based on the parent's true PDG
  int parentPDG = trueParent->PDG;

  // Map PDG codes to category codes
  int categoryCode = 0; // default to "other"

  switch(abs(parentPDG)) {
    case 321:  // K+
      categoryCode = (parentPDG > 0) ? 321 : 321; // Both K+ and K- map to same code
      break;
    case 310:  // K0
      categoryCode = 310;
      break;
    case 211:  // pi+/pi-
      categoryCode = (parentPDG > 0) ? 211 : 211; // Both pi+ and pi- map to same code
      break;
    case 111:  // pi0
      categoryCode = 111;
      break;
    case 2212: // proton
      categoryCode = 2212;
      break;
    case 2112: // neutron
      categoryCode = 2112;
      break;
    case 13:   // mu+/mu-
      categoryCode = (parentPDG > 0) ? 13 : 13; // Both mu+ and mu- map to same code
      break;
    case 11:   // e+/e-
      categoryCode = (parentPDG > 0) ? 11 : 11; // Both e+ and e- map to same code
      break;
    case 22:   // gamma
      categoryCode = 22;
      break;
    default:
      categoryCode = 0; // other
      break;
  }

  anaUtils::_categ->SetCode("vertexparentpdg", categoryCode, CATOTHER);
}

//********************************************************************
bool neutralKaonAnalysis::EventContainsSignalVertices(const AnaEventB& event) const {
//********************************************************************

  // Check if any vertex is classified as signal (k0invtx == 1)
  // We need to check the category system to see if any vertex has been marked as signal

  // Get the number of reconstructed vertex candidates
  const ToyBoxNeutralKaon& neutralKaonBox = static_cast<const ToyBoxNeutralKaon&>(box());

  // If no vertex candidates, no signal
  if (neutralKaonBox.reconVertexCandidates.empty()) {
    return false;
  }

  // Check each vertex candidate to see if it's classified as signal
  // Store signal flags in an array to avoid overwriting
  std::vector<bool> vertexSignalFlags(neutralKaonBox.reconVertexCandidates.size(), false);

  for (size_t vtxIdx = 0; vtxIdx < neutralKaonBox.reconVertexCandidates.size(); vtxIdx++) {
    const auto& vertex = neutralKaonBox.reconVertexCandidates[vtxIdx];
    if (!vertex) continue;

    // Check if this vertex has exactly 2 daughters (K0 -> pi+ pi-)
    if (vertex->NParticles >= 3) { // parent + 2 daughters
      std::vector<AnaTrueParticleB*> daughterTrueParts;

      // Skip the first particle (parent) and collect the daughters
      for (int i = 1; i < vertex->NParticles; i++) {
        AnaParticlePD* recoPart = vertex->Particles[i];
        if (!recoPart) continue;

        // Get the associated true particle
        AnaTrueParticleB* truePart = recoPart->GetTrueParticle();
        if (!truePart) continue;

        daughterTrueParts.push_back(truePart);
      }

      // Check if we have exactly 2 daughters
      if (daughterTrueParts.size() == 2) {
        AnaTrueParticleB* daughter1 = daughterTrueParts[0];
        AnaTrueParticleB* daughter2 = daughterTrueParts[1];

        // Check if both daughters have the same K0 parent (PDG=310 and same ParentID)
        if (daughter1->ParentPDG == 310 && daughter2->ParentPDG == 310) {
          // Check if they have the same parent particle ID (same K0)
          if (daughter1->ParentID == daughter2->ParentID && daughter1->ParentID != 0) {
            // Now check if they are specifically π+ and π- (signal definition)
            if ((daughter1->PDG == 211 && daughter2->PDG == -211) ||
                (daughter1->PDG == -211 && daughter2->PDG == 211)) {

              // Additional condition: Check if the K0 parent has a K+ parent (PDG=321)
              // and that this K+ is the same particle that originated the vertex

              // Find the K0 parent particle in the true particles list
              AnaTrueParticleB* k0Parent = nullptr;
              AnaTrueParticleB** allTrueParticles = event.TrueParticles;
              int nTrueParticles = event.nTrueParticles;

              for (int j = 0; j < nTrueParticles; j++) {
                if (allTrueParticles[j] && allTrueParticles[j]->ID == daughter1->ParentID) {
                  k0Parent = allTrueParticles[j];
                  break;
                }
              }

              if (k0Parent) {
                // Check if the K0 parent has a K+ parent (PDG=321)
                if (k0Parent->ParentPDG == 321) {
                  // Find the K+ parent particle
                  AnaTrueParticleB* kPlusParent = nullptr;
                  for (int k = 0; k < nTrueParticles; k++) {
                    if (allTrueParticles[k] && allTrueParticles[k]->ID == k0Parent->ParentID) {
                      kPlusParent = allTrueParticles[k];
                      break;
                    }
                  }

                  if (kPlusParent) {
                    // Check if this K+ is the same particle that originated the vertex
                    // The vertex parent should be the reconstructed particle corresponding to this K+
                    AnaParticlePD* vertexParent = static_cast<AnaParticlePD*>(vertex->Particles[0]);
                    if (vertexParent) {
                      AnaTrueParticleB* vertexParentTrue = vertexParent->GetTrueParticle();
                      if (vertexParentTrue && vertexParentTrue->ID == kPlusParent->ID) {
                        vertexSignalFlags[vtxIdx] = true; // Mark this vertex as signal
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Check if any vertex is classified as signal
  for (bool isSignal : vertexSignalFlags) {
    if (isSignal) {
      return true; // Found at least one signal vertex
    }
  }

  return false; // No signal vertices found
}
