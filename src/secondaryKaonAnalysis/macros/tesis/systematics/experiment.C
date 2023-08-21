void experiment(){

  Experiment exp("pdsp");

  DataSample* data = new DataSample("/data4/DUNE/migue/analysis/files/data/6-7GeV_prod4a_reco2_microtree_nobranches_2023-04-18.root");
  DataSample* mc   = new DataSample("/data4/DUNE/migue/analysis/files/systematics/6-7GeV_prod4a_microtree_syst_all_0__2023-04-15.root");

  SampleGroup total("total");
  total.AddDataSample(data);
  total.AddMCSample("sys",mc);

  exp.AddSampleGroup("total",total);

  mc->SetCurrentTree("all_syst");

  DrawingTools draw("/data4/DUNE/migue/analysis/files/systematics/6-7GeV_prod4a_microtree_syst_recombination_2023-04-12.root");

  const int nbins = 10;
  double edges[nbins+1] = {0,10,20,40,60,70,80,90,100,150,250};
  draw.Draw(exp,"bestcandidate_chi2_kaon/bestcandidate_chi2_ndf",nbins,edges,"all","accum_level>2","","sys e2");

}
