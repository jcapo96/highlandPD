#include "kaonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  kaonAnalysis* ana = new kaonAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
