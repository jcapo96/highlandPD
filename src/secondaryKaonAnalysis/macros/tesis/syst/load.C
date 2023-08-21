{
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-02-02.root";
  std::string fmc   = dir+"mc/prod4a/6-7GeV_prod4a_microtree_2023-02-02.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* truth = (TTree*)_file1->Get("truth");

  DrawingTools draw(fmc.c_str());
}
