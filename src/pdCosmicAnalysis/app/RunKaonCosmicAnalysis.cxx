#include "kaonCosmicAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  kaonCosmicAnalysis* ana = new kaonCosmicAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
