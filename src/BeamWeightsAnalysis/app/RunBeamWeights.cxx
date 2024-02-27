#include "BeamWeights.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  BeamWeights* ana = new BeamWeights();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
