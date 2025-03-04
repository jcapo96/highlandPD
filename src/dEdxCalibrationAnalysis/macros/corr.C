void corr(){

  TFile* rfile1 = TFile::Open("YZcorrection.root","r");
  TFile* rfile2 = TFile::Open("YZcorrection_alternative_SCE.root","r");

  TH2F* h1 = (TH2F*)rfile1->Get("corr_hist");
  TH2F* h2 = (TH2F*)rfile2->Get("corr_hist");

  TH1F* h_diff = new TH1F("h_diff","h_diff",20,-20,20);
  TH2F* h_diff_map = (TH2F*)h1->Clone();
  h_diff_map->Reset();

  int nbinsz = h_diff_map->GetNbinsX();
  int nbinsy = h_diff_map->GetNbinsY();

  for(int z = 0; z < nbinsz; z++){
    for(int y = 0; y < nbinsy; y++){
      double nom_value = h1->GetBinContent(z+1,y+1);
      double alt_value = h2->GetBinContent(z+1,y+1);
      h_diff->Fill((alt_value-nom_value)/nom_value*100);
      h_diff_map->SetBinContent(z+1,y+1,(alt_value-nom_value)/nom_value*100);
    }
  }
  
  h_diff_map->GetZaxis()->SetRangeUser(0,10);
  h_diff_map->Draw("colz");
}
