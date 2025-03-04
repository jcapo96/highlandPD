#include "secondaryProtonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  secondaryProtonAnalysis* ana = new secondaryProtonAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
