//header
#include "BeamWeights.hxx"

//highland
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"
#include "dEdxUtils.hxx"

//converters
#include "HighlandMiniTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//selections
#include "BeamWeightsSelection.hxx"

//corrections
#include "ParticlePositionSCECorrection.hxx"

//********************************************************************
BeamWeights::BeamWeights(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************
  
}

//********************************************************************
bool BeamWeights::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("BeamWeights.MinAccumLevelToSave"));
                                                       
  // Define standard categories for color drawing
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle

  return true;
}

//********************************************************************
void BeamWeights::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("PDSPAnalyzer",   new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void BeamWeights::DefineSelections(){
//********************************************************************

  sel().AddSelection("Beam Event Selection", "Beam Event Selection", new BeamWeightsSelection(false));  
}

//********************************************************************
void BeamWeights::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
}

//********************************************************************
void BeamWeights::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void BeamWeights::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void BeamWeights::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  standardPDTree::AddStandardVariables_EventInfo(output());
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
}

//********************************************************************
void BeamWeights::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void BeamWeights::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if(addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  standardPDTree::FillStandardVariables_EventInfo(        output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
  standardPDTree::FillStandardVariables_BeamReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamTrue(         output(), GetSpill().Beam);
}

//********************************************************************
void BeamWeights::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  // if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

}

//********************************************************************
void BeamWeights::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  // baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void BeamWeights::FillCategories(){
//********************************************************************

  // For the beam 
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart){
    anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
  }

}
