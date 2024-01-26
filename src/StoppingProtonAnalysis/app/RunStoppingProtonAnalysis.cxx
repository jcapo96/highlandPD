#include "StoppingProtonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  StoppingProtonAnalysis* ana = new StoppingProtonAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
