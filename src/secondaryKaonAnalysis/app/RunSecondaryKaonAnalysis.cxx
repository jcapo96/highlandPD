#include "secondaryKaonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  secondaryKaonAnalysis* ana = new secondaryKaonAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
