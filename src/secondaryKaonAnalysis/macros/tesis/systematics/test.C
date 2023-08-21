void test(){

  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_recombination_2023-04-12.root";

  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

    //drawing tools
  DrawingTools draw(fmc.c_str());
  draw.ChangeToProtoDUNEStyle(); //set default protodune style

  TH1F* h = new TH1F("h","h",700,1000,1700);
  for(int i = 0; i < 100; i++){
    std::stringstream ssi;
    ssi << i;
    h->Fill(syst->Draw("bestcandidate_chi2_kaon/bestcandidate_chi2_ndf",("bestcandidate_chi2_kaon>0 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf>0 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf < 1000 && accum_level["+ssi.str()+"][][]>6").c_str(),"goff"));
  }
  
  std::cout << h->GetRMS()/h->GetMean() << std::endl;
  h->Draw();
  
}
