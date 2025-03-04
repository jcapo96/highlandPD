//simple macro to run the coherent fit over a mc microtree.
//display some plots and save the results in a root file

void CoherentFit_MC(){

  CoherentFit* cf_mc = new CoherentFit("/dune/data2/users/miagarc/microtrees/dunesw_v09_56_00d00/6-7GeV_prod4a_microtree_2022-12-30.root");
  cf_mc->CreateCoherentSamples(84);
  
  cf_mc->GenerateTrueMCHistograms(1,31,2,84,2);
  cf_mc->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  cf_mc->SequentialCoherentFit();

  cf_mc->WriteToRootFile("CoherentFit_MC_2022-12-30"); 
}
