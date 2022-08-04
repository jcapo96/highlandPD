#include "cosmicsAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  cosmicsAnalysis* ana = new cosmicsAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
