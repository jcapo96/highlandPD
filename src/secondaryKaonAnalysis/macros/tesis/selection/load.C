{
  std::string hostname = gSystem->HostName();

  std::string dir = "/dune/app/users/miagarc/technical_note/files/";

  std::string fdata = dir+"data/data_dedx.root";
  std::string fmc   = dir+"mc/mc_dedx.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  // TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

  DrawingTools draw(fmc.c_str());
}
