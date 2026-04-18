// Skim-only: produces pionMomentumSkim.MlSkimOutputFile (tree mlskim). No main analysis microTree.
// Example:
//   pionMomentumSkim.exe -r -o /tmp/unused.root -p .../pionMomentumSkim.parameters.dat mc.txt
// (-o is required by AnalysisLoop; -r prevents creating/filling that file.)

#include "pionMomentumSkimAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char* argv[]) {
  pionMomentumSkimAnalysis* ana = new pionMomentumSkimAnalysis();
  AnalysisLoop loop(ana, argc, argv);
  loop.Execute();
}
