{

  
  const Int_t nbins_thx=5;
  const Int_t nbins_thy=10;
  Float_t thx_min=50;
  Float_t thx_max=300;
  Float_t thy_min=-1;
  Float_t thy_max=0;

  Float_t thx_binsize= (thx_max-thx_min)/(1.*nbins_thx);
  Float_t thy_binsize= (thy_max-thy_min)/(1.*nbins_thy);


  
  TH1F* h[nbins_thx];

  TH1F hh0("hh0","hh0",1,-1,1);
  TH1F hh("hh","hh",1,-1,1);


  TH1F eff("eff","eff",1,thy_min,thy_max);
  eff.GetYaxis()->SetRangeUser(0,0.5);

  eff.Draw();

  c1->Update();
  gPad->Update();

  c1->WaitPrimitive();

  
  for (int i=0;i<nbins_thx;i++){

    stringstream si;
    si << i;
    std::string ssi="h"+si.str();

    
    h[i] = new TH1F(ssi.c_str(),ssi.c_str(),nbins_thy,thy_min,thy_max);

    for (int j=0;j<nbins_thy;j++){

      Float_t thx=thx_min + thx_binsize*(i+0.5);
      Float_t thy=thy_min + thy_binsize*(j+0.5);

      
      std::cout << thx << " " << thy << std::endl;
      
      stringstream ssx;
      ssx << thx;
      string cut_thx = "("+ssx.str()+")";

      stringstream ssy;
      ssy << thy;
      string cut_thy = "("+ssy.str()+")";


      //      std::string cut1 = "abs(cosmics_dir[][0]-"+ cut_thx+")<0.2";

      std::string bcut= "cosmics_length>100 && abs(cosmics_dir[][0])<0.1"; 
      
      std::string cut1 = "abs(cosmics_pos[][0]-"+ cut_thx+")<25";
      std::string cut2 = "abs(cosmics_dir[][1]-"+ cut_thy+")<0.2";
      
      std::string nobeam="(abs(cosmics_pos[][1]-425)>40 || cosmics_pos[][2]>100) && cosmics_pos[][1]<550 && cosmics_pos[][2]<222";

      
      std::string cut0 = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]<236";
      std::string cut  = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]<235 && cosmics_endpos[][2]>232";
      
      //    if (bc>0){
      cut0 = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]>222 ";
      cut  = nobeam+" && "+ bcut + " && "+cut1 + " && " + cut2 + " && cosmics_endpos[][2]>222 && cosmics_endpos[][2]<234";
      //    }

      std::cout << cut0 << std::endl;
      /*      
      draw.Draw(ana,"0",1,-1,1,"all",cut, "","");
      Float_t n  = 1.*draw.GetLastHisto()->GetEntries();
      draw.Draw(ana,"0",1,-1,1,"all",cut0,"","");
      Float_t n0 = 1.*draw.GetLastHisto()->GetEntries();
      */
      ana->Project("hh0","0",cut0.c_str(),"");
      ana->Project("hh","0",cut.c_str(),"");
      Float_t n0 = hh0.GetEntries();
      Float_t n  =  hh.GetEntries();
      
      Float_t w = n/n0;
      Float_t w_e = sqrt(n/(n0*n0)+n*n/(n0*n0*n0));
      
      h[i]->SetBinContent(j+1,w);
      h[i]->SetBinError(j+1,w_e);        
      
      std::cout << n0 << " " << n << " " << w << " +- " << w_e << std::endl;
      
    }
    
    h[i]->SetXTitle("cos #theta");
    h[i]->SetYTitle("track breaking probability");

    h[i]->SetLineColor(i+1);
    
    if (i==0)
      h[i]->Draw("e1 same");
    else
      h[i]->Draw("e1 same");

    c1->WaitPrimitive();
  }
}
