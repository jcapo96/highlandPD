//simple macro to run the coherent fit over a mc microtree.
//display some plots and save the results in a root file

void CoherentFit_MC(){

  gStyle->SetOptStat(0);
  
  CoherentFit* cf_mc = new CoherentFit("../../../../test.root");
  cf_mc->CreateCoherentSamples(84);
  
  cf_mc->GenerateTrueMCHistograms(1,31,2,84,2);
  cf_mc->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  cf_mc->SequentialCoherentFit();

}
