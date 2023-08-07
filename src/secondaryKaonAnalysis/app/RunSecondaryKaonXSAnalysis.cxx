#include "secondaryKaonXSAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  secondaryKaonXSAnalysis* ana = new secondaryKaonXSAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
