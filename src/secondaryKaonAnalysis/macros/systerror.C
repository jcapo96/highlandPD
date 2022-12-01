void systerror(){

  //open root file
  TFile* f = TFile::Open("syst_histos_dQdx_YZcal.root");
  TH2F* h = (TH2F*)f->Get("h_toy_dEdx_RR");

  TH1D* dummy = new TH1D("dummy","dummy",2000,0,20);

  std::vector<double> rel,rr;
  rel.clear();
  rr.clear();

  for(int ibin = 0; ibin < h->GetNbinsX(); ibin++){
    dummy = h->ProjectionY("dummy",ibin+1,ibin+2);
    rel.push_back(dummy->GetRMS()/dummy->GetMean()*100);
    rr.push_back(dummy->GetBinCenter(ibin+1));
  }

  TGraph* tg = new TGraph(rel.size(),&rr[0],&rel[0]);
  tg->GetYaxis()->SetTitle("Relative Error [%]");
  tg->GetXaxis()->SetTitle("Residual range [cm]");
  tg->SetLineColor(2);
  tg->SetLineWidth(3);
  tg->Draw("al");
}

  
