#include "StoppingProtonAnalysis.hxx"
#include "Parameters.hxx"
#include "StoppingProtonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"

#include "pdAnalysisUtils.hxx"

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "kaonTree.hxx"

#include "ParticlePositionSCECorrection.hxx"

#include "baseToyMaker.hxx"

//********************************************************************
StoppingProtonAnalysis::StoppingProtonAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version
  //  ND::versioning().AddPackage("StoppingProtonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
}

//********************************************************************
bool StoppingProtonAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */

  // Initialize the baseAnalysis
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("StoppingProtonAnalysis.MinAccumLevelToSave"));

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("bestcandidate");

  return true;
}

//********************************************************************
void StoppingProtonAnalysis::DefineInputConverters(){
//********************************************************************
  
  /* In this method we add the to the InputManager (accessed by input() ) the InputConverters created 
     in separet files (see, for example, pdIO/PDSPAnalyzerTreeConverter.hxx)
     which define the allowed input file formats.
  */
  
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void StoppingProtonAnalysis::DefineSelections(){
//********************************************************************

  sel().AddSelection("StoppingProtonSelection", "Stopping Proton Selection", new StoppingProtonSelection(false)); // true/false for forcing break  
}

//********************************************************************
void StoppingProtonAnalysis::DefineCorrections(){
//********************************************************************

  // Some corrections are defined in baseAnalysis
  baseAnalysis::DefineCorrections();
  corr().AddCorrection(0, "sce geometric correction", new ParticlePositionSCECorrection());
}

//********************************************************************
void StoppingProtonAnalysis::DefineSystematics(){
//********************************************************************
  
  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void StoppingProtonAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void StoppingProtonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis (these methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  // beam instrumentation info
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  
  // candidate track (beam track) info
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());

  //bestcandidate variables (for comparative with kaons)
  kaonTree::AddKaonVariables_KaonBestCandidateReco    (output());
  kaonTree::AddKaonVariables_KaonBestCandidateTrue    (output());
}

//********************************************************************
void StoppingProtonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void StoppingProtonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);

  // ---------- Additional candidate variables --------------
  if(box().MainTrack){        
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack, beamPart);
    standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);    
    kaonTree::FillKaonVariables_KaonBestCandidateReco    (output(), box().MainTrack);
    kaonTree::FillKaonVariables_KaonBestCandidateTrue    (output(), box().MainTrack);
  }
}

//********************************************************************
void StoppingProtonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

   // Fill the common variables
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase); 
}

//********************************************************************
bool StoppingProtonAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways for the moment
  return true;
}

//********************************************************************
void StoppingProtonAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables
  baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void StoppingProtonAnalysis::FillCategories(){
//********************************************************************

  // For the candidate
  if(box().MainTrack){
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,"");
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,"bestcandidate");
  }
  
  // For the beam track
  AnaParticleB* beam = static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle;
  if(beam)anaUtils::FillCategories(&GetEvent(), beam,"beam");
}
