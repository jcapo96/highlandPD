#include "QPIXAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  QPIXAnalysis* ana = new QPIXAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
