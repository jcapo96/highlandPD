#include "BrokenTracks.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  BrokenTracks* ana = new BrokenTracks();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
