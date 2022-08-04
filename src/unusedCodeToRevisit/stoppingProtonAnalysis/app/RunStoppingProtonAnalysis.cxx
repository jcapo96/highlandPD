#include "stoppingProtonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  stoppingProtonAnalysis* ana = new stoppingProtonAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
