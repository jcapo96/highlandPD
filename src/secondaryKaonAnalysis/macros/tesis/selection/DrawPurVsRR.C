TGraphAsymmErrors* GetPurity(TTree* t){
  
  std::vector<double> eff,eff_l,eff_h,x;

  std::cout << "PURITY COUNTING" << std::endl;

  TH1F* dummy = new TH1F("dummy","dummy",1,0,1);
  for(int i = -1; i < 7; i++){
    double num,den;
    std::stringstream ssi;
    ssi << i;

    std::string cut_den = "accum_level[0][]>"+ssi.str();
    t->Project("dummy","0.5",cut_den.c_str());
    den = dummy->Integral();

    dummy->Reset();

    std::string cut_num = "candidates_truepdg[]==321 && accum_level[0][]>"+ssi.str();
    t->Project("dummy","0.5",cut_num.c_str());
    num = dummy->Integral();

    dummy->Reset();

    eff.push_back(num/den);
    eff_l.push_back(eff.back()-TEfficiency::ClopperPearson(den,num,0.683,false));
    eff_h.push_back(TEfficiency::ClopperPearson(den,num,0.683,true)-eff.back());
    x.push_back(i+1.5);
    std::cout << num << " " << den << " " << eff.back() << std::endl;
  }

  return new TGraphAsymmErrors(x.size(),&x[0],&eff[0],0,0,&eff_l[0],&eff_h[0]);
}

TGraphAsymmErrors* GetEfficiency(TTree* t){
  
  std::vector<double> eff,eff_l,eff_h,x;

  //compute denominator
  TH1F* dummy = new TH1F("dummy","dummy",1,0,1);
  t->Project("dummy","0.5","accum_level[0][0]>-1 && truekaon_truepdg==321");
  double den = dummy->Integral();

  std::cout << "EFFICIENCY COUNTING" << std::endl;
  
  for(int i = -1; i < 7; i++){
    dummy->Reset();
    std::stringstream ssi;
    ssi << i;

    std::string cut_num = "truekaon_truepdg==321 && accum_level[0][0]>"+ssi.str();
    t->Project("dummy","0.5",cut_num.c_str());
    double num = dummy->Integral();

    eff.push_back(num/den);
    eff_l.push_back(eff.back()-TEfficiency::ClopperPearson(den,num,0.683,false));
    eff_h.push_back(TEfficiency::ClopperPearson(den,num,0.683,true)-eff.back());
    x.push_back(i+1.5);
    std::cout << num << " " << den << " " << eff.back() << std::endl;
  }

  return new TGraphAsymmErrors(x.size(),&x[0],&eff[0],0,0,&eff_l[0],&eff_h[0]);
}


void DrawPurVsRR(){

  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  
  //get trees
  std::string dir = "/dune/app/users/miagarc/technical_note/files/";
  std::string fmc   = dir+"mc/mc_dedx.root";

  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* mc = (TTree*)_file1->Get("ana");
  
  gStyle->SetPadRightMargin(0.15);
  TCanvas* c = new TCanvas("c","c",800,500);

  int codes[10]       = {13,11,-211,-321,-13,-11,211,321,2212,999};
  int color[10]       = {2 ,3 ,4   ,94  ,7  ,6  ,31 ,92 ,8   ,48 };
  std::string leg[10] = {"#mu^{-}", "e^{-}", "#pi^{-}", "k^{-}", "#mu^{+}", "e^{+}", "#pi^{+}", "k^{+}" , "p"  , "other"};

  TH1F* h_tot = new TH1F("h_tot","h_tot",30,0,300);
  mc->Project("h_tot","bestcandidate_hit_resrange[]","accum_level[0]>2");
  int nt = h_tot->GetEntries();

  THStack* hs = new THStack("hs","");
  TH1F* h[10];
  int na = 0;
  for(int i = 0; i < 10; i++){
    std::stringstream ssi, ssc;
    ssi << i;
    ssc << codes[i];
    h[i] = new TH1F(("h"+ssi.str()+"").c_str(),("h"+ssi.str()+"").c_str(),30,0,300);

    mc->Project(h[i]->GetName(),"bestcandidate_hit_resrange[]",("accum_level[0]>2 && (bestcandidateparticle=="+ssc.str()+")").c_str());
    na += h[i]->GetEntries();
    h[i]->Divide(h_tot);
    h[i]->Scale(100);
    h[i]->SetFillColor(color[i]);
    h[i]->SetLineColor(color[i]);
    h[i]->SetMarkerColor(color[i]);
    h[i]->SetMarkerSize(0);
    h[i]->SetTitle(leg[i].c_str());
    h[i]->SetName(leg[i].c_str());
    hs->Add(h[i],"histo");
  }

  // h[0]->GetXaxis()->SetTitle("Residual Range [cm]");
  // h[0]->GetYaxis()->SetTitle("Normalized hit purity (%)");
  // h[0]->GetXaxis()->CenterTitle();
  // h[0]->GetYaxis()->CenterTitle();

  hs->Draw("histo");
  hs->GetXaxis()->SetTitle("Residual Range [cm]");
  hs->GetYaxis()->SetTitle("Normalized hit purity (%)");
  hs->GetXaxis()->CenterTitle();
  hs->GetYaxis()->CenterTitle();
  hs->Draw();

  gPad->BuildLegend(0.89,0.2,0.98,0.8);

  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
}


// Types for 'bestcandidateparticle' category
// #mu^{-}          code 13         color 2
// e^{-}          code 11         color 3
// #pi^{-}          code -211         color 4
// k^{-}          code -321         color 94
// #mu^{+}          code -13         color 7
// e^{+}          code -11         color 6
// #pi^{+}          code 211         color 31
// k^{+}          code 321         color 92
//                         p          code 2212         color 8
//                     other          code 999         color 48
//                  no truth          code -1         color 92
//                  sand #mu          code 777         color 51
// -------------------------------------------
