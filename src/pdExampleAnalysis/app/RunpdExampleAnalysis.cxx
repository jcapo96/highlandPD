#include "pdExampleAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  pdExampleAnalysis* ana = new pdExampleAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
