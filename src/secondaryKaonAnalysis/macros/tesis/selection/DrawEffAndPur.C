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


void DrawEffAndPur(){

  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.14);

  TCanvas* c = new TCanvas("c","c",900,500);
  
  //get trees
  std::string dir = "/dune/app/users/miagarc/technical_note/files/";
  std::string fmc   = dir+"mc/mc.root";

  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* truth = (TTree*)_file1->Get("truth");
  
  TGraphAsymmErrors* Pur = GetPurity(mc);

  TGraphAsymmErrors* Eff = GetEfficiency(truth);

  Pur->SetMarkerStyle(20);
  Pur->SetMarkerColor(9);
  Eff->SetMarkerStyle(20);
  Eff->SetMarkerColor(98);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(Pur,"lp");
  mg->Add(Eff,"lp");

  TH1F* h = new TH1F("","",8,0,8);
  h->GetYaxis()->SetRangeUser(0.001,1.5);

  std::string axis[] = {"NO CUT","BEAM TRACK EXISTS","CANDIDATE EXISTS","DAUGHTER TRACK-LIKE","DAUGHTER #chi^{2}","DAUGHTER MOM","ANGLE WRT DAUGHTER","DISTANCE TO DAUGHTER"};
  for(int i = 0; i < 8; i++)h->GetXaxis()->SetBinLabel(i+1,axis[i].c_str());
  h->GetXaxis()->SetTickLength(0);

  h->Draw();

  mg->Draw("p");

  TLine* l[2][8];

  for(int i = 0; i < 8; i++){
    l[0][i] = new TLine(0.15+0.1*i,0.17,0.15+0.1*i,0.19);
    l[1][i] = new TLine(0.15+0.1*i,0.905,0.15+0.1*i,0.925);

    l[0][i]->SetNDC();
    l[1][i]->SetNDC();
	   
    l[0][i]->Draw();
    l[1][i]->Draw();
  }

  //draw legend now                                                                                                                                                                                        
  TLegend* lg = new TLegend(0.08,0.40,0.28,0.60);
  lg->AddEntry(Pur,"Purity","pl");
  lg->AddEntry(Eff,"Efficiency","pl");
  lg->Draw("same");

  gPad->RedrawAxis();


  gPad->SetLogy();


  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.06,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->Update();

  gPad->Print("plots/selection_pureff.pdf");
}
