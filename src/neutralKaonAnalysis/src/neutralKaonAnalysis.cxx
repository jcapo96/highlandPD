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
  standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),20);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),20);
  standardPDTree::AddStandardVariables_BeamTruthDaughters(output(),20);


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
  // Function in standardPDTree.cxx where the truth tree variables are defined: momentum, pdg, etc.
  // Function in standardPDTree.cxx -> beamParticleTruthDaughters()
}

//********************************************************************
void neutralKaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from pdBaseAnalysis (run, event, ...)
  if (addBase) pdBaseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis
  // standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  // standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);

  // ---------- Additional candidate variables --------------
  if(box().MainTrack){
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    standardPDTree::FillStandardVariables_EventInfo(                  output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
    standardPDTree::FillStandardVariables_BeamInstrumentationReco(    output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamInstrumentationTrue(    output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamParticleReco(           output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleTrue(           output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleHitsReco(       output(), box().MainTrack);
    int ndau = std::min(20,(int)box().MainTrack->Daughters.size());
    for(int i = 0; i < ndau; i++){
      standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
      standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
      output().IncrementCounter(standardPDTree::seltrk_ndau);
    }

    AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
    if(trueBeamPart){
      // std::cout << "trueBeamPart->Daughters.size() = " << trueBeamPart->Daughters.size() << std::endl;
      int ndau_truth = std::min(20,(int)trueBeamPart->Daughters.size());
      for(int i = 0; i < ndau_truth; i++){
        AnaTrueParticlePD* truthdau = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles, trueBeamPart->Daughters[i]);
        // int ndaudau_truth = std::min(20,(int)truthdau->Daughters.size());
        // for(int j = 0; j < ndaudau_truth; j++){
        //   AnaTrueParticlePD* truthdau_dau = pdAnaUtils::GetTrueParticle(    GetSpill().TrueParticles
        //   , truthdau->Daughters[j]);
        //   standardPDTree::FillStandardVariables_BeamTruthDaughters(     output(), static_cast<AnaTrueParticlePD*>(truthdau_dau));
        // }
        standardPDTree::FillStandardVariables_BeamTruthDaughters(     output(), static_cast<AnaTrueParticlePD*>(truthdau));
        output().IncrementCounter(standardPDTree::seltrk_truthdau_ndau);
      }
  }
    // standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack, beamPart);
    // standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);
    //standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);
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
  if (part->PDG != 310) return false;
  else return true;
}

//********************************************************************
void neutralKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************
    // Fill the common variables
    pdBaseAnalysis::FillTruthTree(part);
    neutralKaonTree::FillNeutralKaonVariables_TrueNeutralKaonCandidates(output(), &part);
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
