#include "pdCosmicAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  pdCosmicAnalysis* ana = new pdCosmicAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
