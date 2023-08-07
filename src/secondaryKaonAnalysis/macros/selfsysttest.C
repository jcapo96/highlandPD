void selfsysttest(){

  gStyle->SetOptStat(0);

  CoherentFit* cf_mc = new CoherentFit("/home/migue/Documents/PhD/tesis/files/microtrees/systematics/6-7GeV_prod4a_microtree_syst_beamMom_4-00_1toy_2023-04-06.root",true);
  cf_mc->CreateCoherentSamples(999);
  cf_mc->GenerateTrueMCHistograms(1,61,2,999,1,20);
  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* ts_sample = cf_mc->GetTrueSignalSample();
  CoherentSample* tb_sample = cf_mc->GetTrueBackgroundSample();

  ts_sample->SequentialCoherentFit();
  tb_sample->SequentialCoherentFit();
  sb_sample->GetInitialParValuesForCoherentFit();
  
  cf_mc->ComputeSelfSystematicError();
}

