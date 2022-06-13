void systtest(){

  gStyle->SetOptStat(0);
  
  CoherentFit* cf_mc = new CoherentFit("/home/miguel/Public/highland/files/mc/systematics/6GeV-7GeV_prod4a_filtered_syst_2022-05-30.root");
  cf_mc->CreateCoherentSamples(84);
  
  cf_mc->GenerateTrueMCHistograms(1,5,2,84,2);
  cf_mc->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  cf_mc->SequentialCoherentFit();
  cf_mc->PropagateSystematicErrors(true,true);
}

