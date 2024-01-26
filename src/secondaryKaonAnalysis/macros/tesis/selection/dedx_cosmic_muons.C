void dedx_cosmic_muons(){

  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  gStyle->SetOptStat(0);

  TFile* dfile = TFile::Open("/dune/app/users/miagarc/technical_note/files/data/data_dedx.root");
  TTree* d = (TTree*)dfile->Get("ana");

  TH1F* h1 = new TH1F("h1","h1",200,0,20);

  d->Project("h1","bestcandidate_dau_hit_dedx","bestcandidate_chi2_kaon>0 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50");

  double max1 = h1->GetMaximum();
  for(int i = 0; i < h1->GetNbinsX(); i++){
    h1->SetBinContent(i+1,h1->GetBinContent(i+1)/max1);
    h1->SetBinError(i+1,h1->GetBinError(i+1)/max1);
  }
  h1->SetMarkerColor(9);
  h1->SetLineColor(9);

  TFile* cfile = TFile::Open("/dune/app/users/miagarc/larsoft/dunesw_v09_63_01d00_test/calibration/6GeV/5770/test/test.root");
  TH1F* h2 = (TH1F*)cfile->Get("cosmic_muons");

  double max2 = h2->GetMaximum();
  for(int i = 0; i < h2->GetNbinsX(); i++){
    h2->SetBinContent(i+1,h2->GetBinContent(i+1)/max2);
    h2->SetBinError(i+1,h2->GetBinError(i+1)/max2);
  }
  h2->SetMarkerSize(0);
  h2->SetMarkerColor(98);
  h2->SetLineColor(98);

  h1->SetName("");
  h2->SetName("");
  h1->SetTitle("");
  h2->SetTitle("");

  h2->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h2->GetYaxis()->SetTitle("Relative frequency");
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();

  h2->Draw("histo");
  h1->Draw("samee");

  h1->SetTitle("Candidate's Daughter");
  h2->SetTitle("Cosmic rays < 60 cm of RR, run 5770");
  gPad->BuildLegend(0.3,0.5,0.9,0.9,"","pl");

  gPad->RedrawAxis();
  gPad->Update();
  
}

// KE: [
//   10,  
//   14,  
//   20,  
//   30,  
//   40,  
//   80,  
//   100, 
//   140, 
//   200, 
//   300, 
//   400, 
//   800, 
//   1000
// ]

// Range: [
// 0.70437,
// 1.27937,
// 2.37894,
// 4.72636,
// 7.5788,
// 22.0917,
//  30.4441,
//  48.2235,
//  76.1461,
//  123.567,
//  170.845,
//  353.438,
//   441.476
// ]
