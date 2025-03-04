#include "dQdxYZCalibration.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  dQdxYZCalibration* ana = new dQdxYZCalibration();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
