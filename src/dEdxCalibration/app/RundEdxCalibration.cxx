#include "dEdxCalibration.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  dEdxCalibration* ana = new dEdxCalibration();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
