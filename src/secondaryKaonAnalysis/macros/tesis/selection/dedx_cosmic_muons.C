void dedx_cosmic_muons(){

  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  gStyle->SetOptStat(0);

  TFile* dfile = TFile::Open("/exp/dune/app/users/miagarc/technical_note/files/data/data_dedx.root");
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

  TFile* cfile = TFile::Open("/exp/dune/app/users/miagarc/larsoft/dunesw_v09_63_01d00_test/calibration/6GeV/5770/test/test.root");
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
  h2->SetTitle("Cosmic rays < 55 cm of RR, run 5770");
  gPad->BuildLegend(0.3,0.5,0.9,0.9,"","pl");

  gPad->RedrawAxis();
  gPad->Update();
 
  gStyle->SetOptTitle(0);
  TCanvas* c2 = new TCanvas("c2","c2",700,700);
  h1->GetXaxis()->SetRangeUser(1,10);
  h2->GetXaxis()->SetRangeUser(1,10);
  TRatioPlot* r = new TRatioPlot(h2,h1,"diffsig");
  r->Draw();
  r->GetLowerRefYaxis()->SetRangeUser(-20,10);
  r->GetLowerRefYaxis()->SetTitle("(Diff)/error_{Diff}");
  r->GetLowerRefYaxis()->CenterTitle();
  r->GetLowerRefXaxis()->CenterTitle();
  r->GetLowerRefYaxis()->SetTitleSize(r->GetLowerRefYaxis()->GetTitleSize()*5./7);
  r->GetLowerRefXaxis()->SetTitleSize(r->GetLowerRefXaxis()->GetTitleSize()*5./7);
  r->GetUpperRefYaxis()->SetTitleSize(r->GetUpperRefYaxis()->GetTitleSize()*5./7);
  r->GetLowerRefYaxis()->SetLabelSize(r->GetLowerRefYaxis()->GetLabelSize()*5./7);
  r->GetLowerRefXaxis()->SetLabelSize(r->GetLowerRefXaxis()->GetLabelSize()*5./7);
  r->GetUpperRefYaxis()->SetLabelSize(r->GetUpperRefYaxis()->GetLabelSize()*5./7);
  r->GetLowerRefYaxis()->SetTitleOffset(r->GetLowerRefYaxis()->GetTitleOffset()*7./5);
  r->GetLowerRefXaxis()->SetTitleOffset(r->GetLowerRefXaxis()->GetTitleOffset()*6./5);
  r->GetUpperRefYaxis()->SetTitleOffset(r->GetUpperRefYaxis()->GetTitleOffset()*7./5);
  r->GetLowerRefYaxis()->SetLabelOffset(r->GetLowerRefYaxis()->GetLabelOffset()*5./7);
  r->GetLowerRefXaxis()->SetLabelOffset(r->GetLowerRefXaxis()->GetLabelOffset()*5./7);
  r->GetUpperRefYaxis()->SetLabelOffset(r->GetUpperRefYaxis()->GetLabelOffset()*5./7);

  r->GetUpperPad()->cd();

  gPad->Update();
  gPad->BuildLegend(0.3,0.5,0.85,0.85);
  gPad->RedrawAxis();
 
  TLatex tt1;
  tt1.SetNDC();
  TLatex tt2;
  tt2.SetNDC();
  tt2.SetTextAlign(31);
  tt1.DrawLatex(0.10,0.905,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.905,"Data");

  c2->cd();
  gPad->Print("plots/daughter_cosmics_dedx.pdf");
}
