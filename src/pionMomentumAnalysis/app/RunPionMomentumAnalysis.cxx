#include "pionMomentumAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char* argv[]) {
  pionMomentumAnalysis* ana = new pionMomentumAnalysis();
  AnalysisLoop loop(ana, argc, argv);
  loop.Execute();
}
