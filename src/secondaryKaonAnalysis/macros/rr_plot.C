void rr_plot(){

  TFile* file = TFile::Open("../../../test.root");
  if(!file)exit(1);

  TTree* t = (TTree*)file->Get("all_syst");
  if(!t)exit(1);

  TH1F* h0 = new TH1F("h0","h0",40,0,2);
  TH1F* h1 = new TH1F("h1","h1",40,0,2);
  TH1F* h2 = new TH1F("h2","h2",40,0,2);
  TH1F* h3 = new TH1F("h3","h3",40,0,2);
  
  t->Project("h0","bestcandidate_hit_resrange[0]","bestcandidate_hit_resrange[0]>-1 && accum_level[0][]>8");
  t->Project("h1","bestcandidate_hit_resrange_toy[0][0]","bestcandidate_hit_resrange_toy[0][0]>-1 && accum_level[0][]>8");
  t->Project("h2","bestcandidate_hit_resrange_toy[1][0]","bestcandidate_hit_resrange_toy[1][0]>-1 && accum_level[1][]>8");
  t->Project("h3","bestcandidate_hit_resrange_toy[2][0]","bestcandidate_hit_resrange_toy[2][0]>-1 && accum_level[2][]>8");
  
  h0->SetLineColor(1);
  h1->SetLineColor(2);
  h2->SetLineColor(3);
  h3->SetLineColor(4);

  h0->SetLineWidth(3);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
  h3->SetLineWidth(3);

  h0->Draw();
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
}


