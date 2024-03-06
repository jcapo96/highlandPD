#include "dQdxXCalibration.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  dQdxXCalibration* ana = new dQdxXCalibration();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
