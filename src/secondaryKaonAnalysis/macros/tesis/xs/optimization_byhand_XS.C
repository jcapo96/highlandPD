//open tree as global to be called anywhere
TFile *rfile = TFile::Open("/data4/DUNE/migue/analysis/files/mc/prod4a/xs_corrected_03-05-2023.root");
TTree* mc  = (TTree*)rfile->Get("ana");

//define cuts globally
std::string signal_cut   = "(seltrk_dau_truepdg==321)";
std::string chi2_cut   = "(seltrk_dau_chi2_kaon>0 && seltrk_dau_chi2_kaon/seltrk_dau_chi2_ndf<3.65)";
std::string pre_cut = "seltrk_dau_chi2_kaon>0";//chi2_cut;

//*********************************************************
double GetSignalBefore(){
//*********************************************************

  TH1F* h1 = new TH1F("h1","h1",1,0,1);

  mc->Project("h1","0.5",(""+signal_cut+" && "+pre_cut+"").c_str());

  return h1->GetEntries();

  // return mc->Draw("0.5",
  // 		  (""+signal_cut+" && "+pre_cut+"").c_str(),
  // 		  "goff"); 
}

//*********************************************************
double GetSignal(double a, double b){
//*********************************************************
  std::stringstream ssa, ssb;
  ssa << a;
  ssb << b;
  std::string lower_cut = "seltrk_dau_chi2_kaon/seltrk_dau_chi2_ndf>"+ssa.str();
  std::string upper_cut = "seltrk_dau_chi2_kaon/seltrk_dau_chi2_ndf<"+ssb.str();

  TH1F* h2 = new TH1F("h2","h2",1,0,1);

  mc->Project("h2","0.5",(""+signal_cut+" && "+pre_cut+" && "+lower_cut+" && "+upper_cut+"").c_str());

  return h2->GetEntries();
  // return mc->Draw("0.5",
  // 		  (""+signal_cut+" && "+pre_cut+" && "+lower_cut+" && "+upper_cut+"").c_str(),
  // 		  "goff");
}

//*********************************************************
double GetAll(double a, double b){
//*********************************************************
  std::stringstream ssa, ssb;
  ssa << a;
  ssb << b;
  std::string lower_cut = "seltrk_dau_chi2_kaon/seltrk_dau_chi2_ndf>"+ssa.str();
  std::string upper_cut = "seltrk_dau_chi2_kaon/seltrk_dau_chi2_ndf<"+ssb.str();

  TH1F* h3 = new TH1F("h3","h3",1,0,1);

  mc->Project("h3","0.5",(""+pre_cut+" && "+lower_cut+" && "+upper_cut+"").c_str());

  return h3->GetEntries();

  // return mc->Draw("0.5",
  // 		  (""+pre_cut+" && "+lower_cut+" && "+upper_cut+"").c_str(),
  // 		  "goff");
}

//*********************************************************
double GetEfficiency(double a, double b){
//*********************************************************

  double signal_before = GetSignalBefore();
  double signal_after  = GetSignal(a,b);
  return signal_after/signal_before;
}

//*********************************************************
double GetPurity(double a, double b){
//*********************************************************

  double signal = GetSignal(a,b);
  double all    = GetAll(a,b);
  return signal/all;
}

//*********************************************************
double EventSelection(double a, double b){
//*********************************************************

  if(a>=b){
    std::cout << "a>b is not allowed!" << std::endl;
    return -9999;
  }
  
  double pur = GetPurity(a,b);
  double eff = GetEfficiency(a,b);

  std::cout << "cut: " << a << "-" << b << " || PUR = " << pur << " || EFF = " << eff << std::endl;
  return pur*eff*10000;
}

//*********************************************************
void optimization_byhand_XS(){
//*********************************************************

  double bmin = 0.5;
  double bmax = 1.0;
  double step = 0.01;
  double b    = bmin;
  
  TGraph* tg = new TGraph();
  int ipoint = 0;

  while(b<bmax){
    std::cout << b << std::endl;
    //tg->SetPoint(ipoint,b,EventSelection(2,b));
    tg->SetPoint(ipoint,b,EventSelection(b,3.65));
    b += step;
    ipoint++;
  }

  tg->SetMarkerStyle(20);
  tg->Draw("ap");
  
}
