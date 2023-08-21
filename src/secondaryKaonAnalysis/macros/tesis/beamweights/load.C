{
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-04-04.root";
  std::string fmc   = dir+"mc/prod4a/6-7GeV_prod4a_microtree_2023-04-04.root";
  //std::string fdata = dir+"data/BeamWeights_6GeV_data.root";
  //std::string fmc   = dir+"mc/prod4a/BeamWeights_6GeV_mc.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  // TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

  DrawingTools draw(fmc.c_str());
}
