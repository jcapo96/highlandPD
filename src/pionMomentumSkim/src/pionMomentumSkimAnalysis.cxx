#include "pionMomentumSkimAnalysis.hxx"
#include "pionMomentumSelection.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
pionMomentumSkimAnalysis::pionMomentumSkimAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************
}

//********************************************************************
pionMomentumSkimAnalysis::~pionMomentumSkimAnalysis() {
//********************************************************************
}

//********************************************************************
bool pionMomentumSkimAnalysis::Initialize() {
//********************************************************************
  if (!baseAnalysis::Initialize()) return false;

  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pionMomentumSkim.MinAccumLevelToSave"));

  if (!_skim.Initialize("pionMomentumSkim")) return false;

  return true;
}

//********************************************************************
void pionMomentumSkimAnalysis::Finalize() {
//********************************************************************
  _skim.Finalize();
}

//********************************************************************
void pionMomentumSkimAnalysis::FinalizeBunch() {
//********************************************************************
  _skim.FillForBunch(GetEvent(), GetSpill(), GetBunch());
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineInputConverters() {
//********************************************************************
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineSelections() {
//********************************************************************
  sel().AddSelection("pionMomentumSelection", "Pion momentum skim (pass all)", new pionMomentumSelection(false));
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineCorrections() {
//********************************************************************
  baseAnalysis::DefineCorrections();
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineSystematics() {
//********************************************************************
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineConfigurations() {
//********************************************************************
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
}

//********************************************************************
void pionMomentumSkimAnalysis::DefineTruthTree() {
//********************************************************************
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void pionMomentumSkimAnalysis::FillMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);
}

//********************************************************************
void pionMomentumSkimAnalysis::FillToyVarsInMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
bool pionMomentumSkimAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* /*part*/) {
//********************************************************************
  return false;
}

//********************************************************************
void pionMomentumSkimAnalysis::FillTruthTree(const AnaTrueParticlePD& /*part*/) {
//********************************************************************
}

//********************************************************************
void pionMomentumSkimAnalysis::FillCategories() {
//********************************************************************
}
