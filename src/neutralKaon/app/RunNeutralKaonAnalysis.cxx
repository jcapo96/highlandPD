#include "neutralKaonAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  neutralKaonAnalysis* ana = new neutralKaonAnalysis();
  AnalysisLoop loop(ana, argc, argv);
  loop.Execute();
}
