#include "dEdxStudies.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  dEdxStudies* ana = new dEdxStudies();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
