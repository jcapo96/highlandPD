#include "neutralKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "neutralKaonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"

#include "pdAnalysisUtils.hxx"
#include "standardPDTree.hxx"

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

// #include "ParticlePositionSCECorrection.hxx"

#include "baseToyMaker.hxx"

//********************************************************************
neutralKaonAnalysis::neutralKaonAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************

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

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("bestcandidate");

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
  // corr().AddCorrection(0, "sce geometric correction", new ParticlePositionSCECorrection());
}

//********************************************************************
void neutralKaonAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in pdBaseAnalysis (highland/src/highland2/pdBaseAnalysis)
  pdBaseAnalysis::DefineSystematics();
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

  // standardPDTree::AddStandardVariables_AllParticlesReco(output(), 10);
  // standardPDTree::AddStandardVariables_AllParticlesTrue(output(), 10);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_EventInfo(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());
  standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
  standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),50);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),50);
  standardPDTree::AddStandardVariables_BeamTruthDaughters(output(),50);

  AddVarF(output(), seltrk_goodk0, "K0 is daughter & K0 daughters are daughters of particle");
  AddVarI(output(), seltrk_dau_trueparentpdg, "Parent PDG of reco daughter");


  // Add standard sets of variables for ProtoDUNE analysis (these methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  // beam instrumentation info
  // standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  // standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());

  // // candidate track (beam track) info
  // standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  // standardPDTree::AddStandardVariables_BeamParticleReco(output());
  //standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
}

//********************************************************************
void neutralKaonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  pdBaseAnalysis::DefineTruthTree();
  neutralKaonTree::AddNeutralKaonVariables_TrueNeutralKaonCandidates(output());
  // Add dynamic daughters variables up to a maximum (only filled for existing daughters)
  std::cout << "[neutralKaonAnalysis] Adding dynamic true daughters (max=200)" << std::endl;
  neutralKaonTree::AddNeutralKaonVariables_TrueDaughtersDynamic(output(), 200);
  std::cout << "[neutralKaonAnalysis] Done adding dynamic true daughters" << std::endl;
  neutralKaonTree::AddNeutralKaonVariables_TrueParentCandidates(output());
  neutralKaonTree::AddNeutralKaonVariables_TrueGrandParentCandidates(output());
  // Function in standardPDTree.cxx where the truth tree variables are defined: momentum, pdg, etc.
  // Function in standardPDTree.cxx -> beamParticleTruthDaughters()
}

//********************************************************************
void neutralKaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  if (addBase) pdBaseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);

  // ---------- Additional candidate variables --------------
  if(box().MainTrack){
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
    standardPDTree::FillStandardVariables_EventInfo(                  output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
    standardPDTree::FillStandardVariables_BeamInstrumentationReco(    output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamInstrumentationTrue(    output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamParticleReco(           output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleTrue(           output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleHitsReco(       output(), box().MainTrack);
    int ndau = (int)box().MainTrack->Daughters.size();
    bool goodk0 = false;
    for(int i = 0; i < ndau; i++){
      // cout << "box().MainTrack->Daughters[i]->TrueObject->ID = " << box().MainTrack->Daughters[i] << endl;
      standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
      standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
      output().FillVar(seltrk_dau_trueparentpdg, static_cast<AnaTrueParticlePD*>(box().MainTrack->Daughters[i]->TrueObject)->ParentPDG);
      if (trueBeamPart) {
        standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack, beamPart);
        standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);
        standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);
        int ndau_truth = (int)trueBeamPart->Daughters.size();
        // cout << "trueBeamPart->Daughters.size() = " << trueBeamPart->Daughters.size() << endl;
        for(int j = 0; j < ndau_truth; j++){
          // cout << "trueBeamPart->Daughters[j] = " << trueBeamPart->Daughters[j] << endl;
          AnaTrueParticlePD* truthdau = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles, trueBeamPart->Daughters[j]);
          if (!truthdau) continue;
          // cout << "truthdau->PDG = " << truthdau->PDG << endl;
          if (truthdau->PDG == 310) {
            int ndau_truthk0 = (int)truthdau->Daughters.size();
            // cout << "truthdau->Daughters.size() = " << truthdau->Daughters.size() << endl;
            // cout << "trueBeamPart->Daughters[j]: " << trueBeamPart->Daughters[j] << endl;
            for(int k = 0; k < ndau_truthk0; k++){
              // cout << "truthdau->Daughters[k]: " << truthdau->Daughters[k] << endl;
              AnaTrueParticlePD* truthdauk0 = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles, truthdau->Daughters[k]);
              // cout << truthdauk0 << endl;
              // cout << "GetSpill.TrueParticles.size(): " << GetSpill().TrueParticles.size() << endl;
              if (!truthdauk0) continue;
              // cout << "truthdauk0->ID: " << truthdauk0->ID << endl;
              // cout << "box().MainTrack->Daughters[i]->TrueObject->ID" << box().MainTrack->Daughters[i]->TrueObject->ID << endl;
              if (truthdauk0->ID == box().MainTrack->Daughters[i]->TrueObject->ID) {
                // cout << "True - K0 has " << ndau_truthk0 << " daughters" << endl;
                goodk0 = true;
              }
              output().IncrementCounter(neutralKaonAnalysis::seltrk_truthk0_ndau);
            }
          }
          output().IncrementCounter(standardPDTree::seltrk_truthdau_ndau);
        }
      }
    float seltrk_goodk0_true = 100.0;
    float seltrk_goodk0_false = 0.0;
    if (goodk0) {
      output().FillVar(seltrk_goodk0, seltrk_goodk0_true);
    }
    else {
      output().FillVar(seltrk_goodk0, seltrk_goodk0_false);
    }
    output().IncrementCounter(standardPDTree::seltrk_ndau);
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
  // Only fill truth tree if there is a particle with PDG=310 in the daughters
  if (!part) return false;

  // Check if any daughter has PDG=310 (neutral kaon)
  const std::vector<int>& daughters = part->Daughters;
  for (size_t i = 0; i < daughters.size(); i++) {
    AnaTrueParticlePD* daughter = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, daughters[i]);
    if (daughter && daughter->PDG == 310) {
      return true; // Found a neutral kaon daughter
    }
  }

  return false; // No neutral kaon daughters found
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************
    // Fill the common variables
    pdBaseAnalysis::FillTruthTree(part);
    neutralKaonTree::FillNeutralKaonVariables_TrueNeutralKaonCandidates(output(), &part);
    AnaTrueParticlePD* truthpar = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles, part.ParentID);
    if (truthpar) {
      neutralKaonTree::FillNeutralKaonVariables_TrueParentCandidates(output(), truthpar);
      AnaTrueParticlePD* truthgpar = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles, truthpar->ParentID);
      if (truthgpar) {
        neutralKaonTree::FillNeutralKaonVariables_TrueGrandParentCandidates(output(), truthgpar);
      }
    }
    if(&part){
      // Fill dynamic daughters (support up to 200, but only fill existing)
      neutralKaonTree::FillNeutralKaonVariables_TrueDaughtersWithCollection(
        output(), &part, GetSpill().TrueParticles, 1000);
    }
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
}
