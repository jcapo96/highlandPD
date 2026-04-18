#include "pionMomentumAnalysis.hxx"
#include "PionMomentumMlSkimInputConverter.hxx"
#include "pionMomentumSelection.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
pionMomentumAnalysis::pionMomentumAnalysis(AnalysisAlgorithm* ana) : pdBaseAnalysis(ana) {
//********************************************************************
}

//********************************************************************
pionMomentumAnalysis::~pionMomentumAnalysis() {
//********************************************************************
}

//********************************************************************
bool pionMomentumAnalysis::Initialize() {
//********************************************************************
  if (!baseAnalysis::Initialize()) return false;

  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pionMomentumAnalysis.MinAccumLevelToSave"));

  return true;
}

//********************************************************************
void pionMomentumAnalysis::Finalize() {
//********************************************************************
}

//********************************************************************
void pionMomentumAnalysis::DefineInputConverters() {
//********************************************************************
  input().AddConverter("mlskim", new PionMomentumMlSkimInputConverter());
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void pionMomentumAnalysis::DefineSelections() {
//********************************************************************
  sel().AddSelection("pionMomentumSelection", "Pion momentum analysis (pass all)", new pionMomentumSelection(false));
}

//********************************************************************
void pionMomentumAnalysis::DefineCorrections() {
//********************************************************************
  baseAnalysis::DefineCorrections();
}

//********************************************************************
void pionMomentumAnalysis::DefineSystematics() {
//********************************************************************
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void pionMomentumAnalysis::DefineConfigurations() {
//********************************************************************
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void pionMomentumAnalysis::DefineMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
}

//********************************************************************
void pionMomentumAnalysis::DefineTruthTree() {
//********************************************************************
  baseAnalysis::DefineTruthTree();
  AddVarI(output(), pion_true_parent_of_secondaries_pdg,
          "True PDG of particle parenting event secondaries (skim seed; expect 211 for pi+ skim)");
}

//********************************************************************
void pionMomentumAnalysis::FillMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);
}

//********************************************************************
void pionMomentumAnalysis::FillToyVarsInMicroTrees(bool addBase) {
//********************************************************************
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
bool pionMomentumAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* /*part*/) {
//********************************************************************
  return input().GetConverter().Name() == "mlskim";
}

//********************************************************************
void pionMomentumAnalysis::FillTruthTree(const AnaTrueParticlePD& part) {
//********************************************************************
  pdBaseAnalysis::FillTruthTree(part);
  if (input().GetConverter().Name() != "mlskim") return;
  output().FillVar(pion_true_parent_of_secondaries_pdg, part.PDG);
}

//********************************************************************
void pionMomentumAnalysis::FillCategories() {
//********************************************************************
}
