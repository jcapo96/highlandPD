{

  const Int_t nbins_thx=10;
  const Int_t nbins_thy=10;
  Float_t thx_min=-1;
  Float_t thx_max=1;
  Float_t thy_min=-1;
  Float_t thy_max=0;

  Float_t thx_binsize= (thx_max-thx_min)/(1.*nbins_thx);
  Float_t thy_binsize= (thy_max-thy_min)/(1.*nbins_thy);


  TH2F* eff0 = new TH2F("eff0","eff0",nbins_thx,thx_min,thx_max,nbins_thy,thy_min,thy_max);
  TH2F* eff1 = new TH2F("eff1","eff1",nbins_thx,thx_min,thx_max,nbins_thy,thy_min,thy_max);


    
  std::string bcut= "cosmics_length>100 && abs(cosmics_dir[][0])<1.1"; 
  
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
  ana->Project("eff0","cosmics_dir[][1]:cosmics_dir[][0]",cut0.c_str());
  ana->Project("eff1","cosmics_dir[][1]:cosmics_dir[][0]",cut.c_str());


  gStyle->SetOptStat(0);
  
  
  TH2F *ratio = (TH2F*)eff1->Clone("ratio"); ratio->Divide(eff1,eff0);
  ratio->GetZaxis()->SetRangeUser(0., 0.3);

  ratio->GetXaxis()->SetTitleSize(0.04);
  ratio->GetYaxis()->SetTitleSize(0.04);
  ratio->SetXTitle("cos #theta_{X}");
  ratio->SetYTitle("cos #theta_{Y}");
  ratio->Draw("colz");
}
