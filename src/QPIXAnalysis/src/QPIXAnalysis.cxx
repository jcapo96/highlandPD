#include "QPIXAnalysis.hxx"
#include "Parameters.hxx"
#include "QPIXSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"

#include "QPIXUtils.hxx"

#include "QPIXTreeConverter.hxx"

#include "baseToyMaker.hxx"

//********************************************************************
QPIXAnalysis::QPIXAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version
  //ND::versioning().AddPackage("QPIXAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("QPIXANALYSISROOT")));
}

//********************************************************************
bool QPIXAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */

  // Initialize the baseAnalysis
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(-1);//ND::params().GetParameterI("QPIXAnalysis.MinAccumLevelToSave"));

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  QPIXUtils::AddCustomCategories();

  return true;
}

//********************************************************************
void QPIXAnalysis::DefineInputConverters(){
//********************************************************************
  
  /* In this method we add the to the InputManager (accessed by input() ) the InputConverters created 
     in separet files (see, for example, pdIO/PDSPAnalyzerTreeConverter.hxx)
     which define the allowed input file formats.
  */
  
  input().AddConverter("minitreefiltered", new QPIXTreeConverter("event_tree"));
}

//********************************************************************
void QPIXAnalysis::DefineSelections(){
//********************************************************************

  sel().AddSelection("QPIXSelection", "Dummy", new QPIXSelection(false)); // true/false for forcing break  
}

//********************************************************************
void QPIXAnalysis::DefineCorrections(){
//********************************************************************
  
  // Some corrections are defined in baseAnalysis
  baseAnalysis::DefineCorrections();
}

//********************************************************************
void QPIXAnalysis::DefineSystematics(){
//********************************************************************
  
  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void QPIXAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void QPIXAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
  QPIXUtils::AddQPIXVariables(output());
}

//********************************************************************
void QPIXAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void QPIXAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  QPIXUtils::FillQPIXVariables(output(),&GetEvent());
}

//********************************************************************
void QPIXAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

   // Fill the common variables
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase); 
}

//********************************************************************
bool QPIXAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways for the moment
  return true;
}

//********************************************************************
void QPIXAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables
  //baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void QPIXAnalysis::FillCategories(){
//********************************************************************

  QPIXUtils::FillCustomCategories(&GetEvent());
}
