const double ZMAX = 240;
const double ZMIN = 0;
const double YMAX = 600;
const double YMIN = 0;
const double STEP = 5;

//********************************************
Double_t Langaus(Double_t *x, Double_t *par) {
//********************************************
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  //Double_t np = 100.0;      // number of convolution steps
  Double_t np = 1000.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

//********************************************
void toy_yz(){
//********************************************
  
  //file name
  std::string filename = "/dune/data/users/miagarc/toy.root";
  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* tree  = (TTree*)rfile->Get("all_syst");

  //histograms and functions for calculations
  TH1F* hsyst = new TH1F("h","h",1000,0,100);
  TH1F* hdummy = new TH1F("hdummy","hdummy",100,0,100);
  TF1* f = new TF1("f",Langaus,0,80,4);
  f->SetParameters(6000,60,5,5);

  double global_mpv, global_stat_error, global_syst_error;
  //get global mean and errors
  std::cout << "global mean calculation" << std::endl;
  for(int itoy = 0; itoy < 100; itoy++){
    std::cout << "toy " << itoy << std::endl;
    std::stringstream sitoy;
    sitoy << itoy;
    tree->Project("hdummy",("toy_hit_dqdx["+sitoy.str()+"][]").c_str(),
		  ("toy_hit_z["+sitoy.str()+"][]>0 && toy_hit_z["+sitoy.str()+"][]<695 && toy_hit_y["+sitoy.str()+"][]>0 && toy_hit_y["+sitoy.str()+"][]<600").c_str(),"");
    //hdummy = (TH1F*)gPad->GetPrimitive("hd");
    hdummy->Fit("f","QNO");
    hsyst->Fill(f->GetParameter(1));
    if(itoy==0)global_stat_error = f->GetParError(1);
    hdummy->Reset();
  }
  hsyst->Draw();
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  double global_rel_syst_error = global_syst_error / global_mpv;
  std::cout << global_mpv << " +/- " << global_stat_error << " +/- " << global_syst_error << std::endl;

  hsyst->Reset();
  int nbinsz = (ZMAX-ZMIN)/STEP;
  int nbinsy = (YMAX-YMIN)/STEP;
  std::vector<double> local_mpv, local_stat_error, local_syst_error;
  local_mpv.clear(); local_stat_error.clear(); local_syst_error.clear();
  TH2F* h2 = new TH2F("h","h",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX);
  for(int iz = 0; iz < nbinsz; iz++){
    std::stringstream szl,szh;
    szl << ZMIN + iz*STEP;
    szh << ZMIN + (iz+1)*STEP;
    for(int iy = 0; iy < nbinsy; iy++){
      std::stringstream syl,syh;
      syl << YMIN + iy*STEP;
      syh << YMIN + (iy+1)*STEP;
      std::cout << "bins " << iz << "/" << nbinsz << " & " << iy << "/" << nbinsz << std::endl;
      for(int itoy = 0; itoy < 100; itoy++){
	std::stringstream sitoy; 
        sitoy << itoy;
	//std::cout << "toy " << itoy << std::endl;

	//tree->Draw(("toy_track_hit_dqdx["+sitoy.str()+"][][]>>hd(100,0,100)").c_str(),
	//	   ("toy_track_hit_z["+sitoy.str()+"][][]>"+szl.str()+" && toy_track_hit_z["+sitoy.str()+"][][]<"+szh.str()+" && toy_track_hit_y["+sitoy.str()+"][][]>"+syl.str()+" && toy_track_hit_y["+sitoy.str()+"][][]<"+syh.str()+"").c_str(),"goff");
	// tree->Project("hdummy",("toy_track_hit_dqdx["+sitoy.str()+"][][]").c_str(),
	// 	   ("toy_track_hit_z["+sitoy.str()+"][][]>"+szl.str()+" && toy_track_hit_z["+sitoy.str()+"][][]<"+szh.str()+" && toy_track_hit_y["+sitoy.str()+"][][]>"+syl.str()+" && toy_track_hit_y["+sitoy.str()+"][][]<"+syh.str()+"").c_str());
	tree->Project("hdummy",("toy_hit_dqdx["+sitoy.str()+"][]").c_str(),
		      ("toy_hit_z["+sitoy.str()+"][]>"+szl.str()+" && toy_hit_z["+sitoy.str()+"][]<"+szh.str()+" && toy_hit_y["+sitoy.str()+"][]>"+syl.str()+" && toy_hit_y["+sitoy.str()+"][]<"+syh.str()+"").c_str());
	//hdummy = (TH1F*)gPad->GetPrimitive("hd");
	hdummy->Fit("f","QNO");
	hsyst->Fill(f->GetParameter(1));
	if(itoy==0)local_stat_error.push_back(f->GetParError(1));
	hdummy->Reset();
      }
      local_mpv.push_back(hsyst->GetMean());
      local_syst_error.push_back(hsyst->GetRMS());
      std::cout << local_mpv.back() << " +/- " << local_stat_error.back() << " +/- " << local_syst_error.back() << std::endl;
      hsyst->Draw();gPad->Update();//gPad->WaitPrimitive();
      hsyst->Reset();
      double local_rel_syst_error = local_syst_error.back()/local_mpv.back();
      h2->SetBinContent(iz+1,iy+1,sqrt(pow(global_rel_syst_error,2)+pow(local_rel_syst_error,2))*100);
      std::cout << sqrt(pow(global_rel_syst_error,2)+pow(local_rel_syst_error,2))*100 << std::endl;
    }
  }
  h2->Draw("colz");

  TFile* wfile = TFile::Open("syst_effect.root","NEW");
  h2->Write();
  wfile->Close();
}
