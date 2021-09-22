{


  TFile *_file0 = TFile::Open("test1.root");
  TFile *_file1 = TFile::Open("6GeV_prod4a_00_01_cosmics.root");

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");


  DrawingTools draw("6GeV_prod4a_00_01_cosmics.root");
  
  const Int_t nbins_thx=10;
  const Int_t nbins_thy=10;
  Float_t thx_min=-1;
  Float_t thx_max=1;
  Float_t thy_min=-1;
  Float_t thy_max=0;

  Float_t thx_binsize= (thx_max-thx_min)/(1.*nbins_thx);
  Float_t thy_binsize= (thy_max-thy_min)/(1.*nbins_thy);


  TH1F* eff0_d  = new TH1F("eff0_d" ,"eff0_d", nbins_thy,thy_min,thy_max);
  TH1F* eff1_d  = new TH1F("eff1_d" ,"eff1_d", nbins_thy,thy_min,thy_max);
  TH1F* eff0_mc = new TH1F("eff0_mc","eff0_mc",nbins_thy,thy_min,thy_max);
  TH1F* eff1_mc = new TH1F("eff1_mc","eff1_mc",nbins_thy,thy_min,thy_max);


    
  std::string bcut= "cosmics_length>100 && abs(cosmics_dir[][0])>0.1 && abs(cosmics_dir[][0])<0.2"; 
  
  std::string cut1 = "(1==1)";//"abs(cosmics_pos[][0]-"+ cut_thx+")<25";
  std::string cut2 = "(1==1)";//"abs(cosmics_dir[][1]-"+ cut_thy+")<0.2";
  
  std::string nobeam="(abs(cosmics_pos[][1]-425)>40 || cosmics_pos[][2]>100) && cosmics_pos[][1]<550 && cosmics_pos[][2]<222";
  
  
  std::string cut0 = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]<236";
  std::string cut  = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]<235 && cosmics_endpos[][2]>232";
  
  //    if (bc>0){
  cut0 = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]>222 ";
  cut  = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]>222 && cosmics_endpos[][2]<234";
  //    }
  
  std::cout << cut0 << std::endl;

  d->Project("eff0_d","cosmics_dir[][1]",cut0.c_str());
  d->Project("eff1_d","cosmics_dir[][1]",cut.c_str());

  mc->Project("eff0_mc","cosmics_dir[][1]",cut0.c_str());
  mc->Project("eff1_mc","cosmics_dir[][1]",cut.c_str());


  eff0_d->Sumw2();
  eff1_d->Sumw2();
  eff0_mc->Sumw2();
  eff1_mc->Sumw2();
    
  gStyle->SetOptStat(0);
  
  TCanvas c("c","c",1400,500);
  c.Divide(2,1);

  c.cd(1);
  //  TH1F *ratio_d = (TH1F*)eff1_d->Clone("ratio_d");
  TGraphAsymmErrors* ratio_d = new TGraphAsymmErrors(eff1_d->GetNbinsX());
  ratio_d->Divide(eff1_d,eff0_d);


  
  //  ratio_d->GetXaxis()->SetTitleSize(0.04);
  ratio_d->GetXaxis()->SetTitle("cos #theta_{Y}");
  ratio_d->GetYaxis()->SetTitle("track breaking probability");
  ratio_d->Draw("AP");


  c.cd(1);
  //  TH1F *ratio_mc = (TH1F*)eff1_mc->Clone("ratio_mc");
  TGraphAsymmErrors* ratio_mc = new TGraphAsymmErrors(eff1_d->GetNbinsX());
  ratio_mc->Divide(eff1_mc,eff0_mc,"cl=0.683 b(1,1) mode");

  //  ratio_mc->GetXaxis()->SetTitleSize(0.04);
  //  ratio_mc->SetXTitle("cos #theta_{Y}");
  ratio_mc->SetLineColor(2);
  ratio_mc->Draw("P same");
}
