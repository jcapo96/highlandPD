const int ZMAX = 250;
const int ZMIN = 0;
const int YMAX = 600;
const int YMIN = 0;
const int STEP = 5;
const int NTOYS   = 100;
const int NHITS   = 1000;

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
bool IsValidHit(const double z,    const double y,
		const double zmin, const double zmax,
		const double ymin, const double ymax){
//********************************************

  return (z > zmin && z < zmax && y > ymin && y < ymax);
}

//********************************************
void toy_yz_branch(){
//********************************************
  
  //file name
  std::string filename = "/dune/data/users/miagarc/toy_lifetime.root";
  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* tree  = (TTree*)rfile->Get("all_syst");

  Float_t hit_y [NTOYS][NHITS]    = {0};
  Float_t hit_z [NTOYS][NHITS]    = {0};
  Float_t hit_dqdx [NTOYS][NHITS] = {0};

  tree->SetBranchAddress("toy_hit_y"   , hit_y   );
  tree->SetBranchAddress("toy_hit_z"   , hit_z   );
  tree->SetBranchAddress("toy_hit_dqdx", hit_dqdx);

  int nentries = tree->GetEntries();
  
  //histograms and functions for calculations
  TH1F* hdummy[NTOYS];
  TH1F* hsyst = new TH1F("hsyst","hsyst",1000,0,100);
  //initialize histograms
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::stringstream stoy;
    stoy << itoy;
    hdummy[itoy] = new TH1F(("h"+stoy.str()+"").c_str(),
			   ("h"+stoy.str()+"").c_str(),
			   100,0,100);
  }

  //get global mean and errors
  double global_mpv, global_stat_error, global_syst_error;
  std::cout << "global mpv calculation" << std::endl;
  //loop over entries
  for(int ientry = 0; ientry < nentries; ientry++){
    tree->GetEntry(ientry);
    if(ientry%20000==0)std::cout << ientry << "/" << nentries << std::endl;
    //loop over toys
    for(int itoy = 0; itoy < NTOYS; itoy++){
      //loop over hits
      for(int ihit = 0; ihit < NHITS; ihit++){
	if(IsValidHit(hit_z[itoy][ihit],hit_y[itoy][ihit],ZMIN,ZMAX,YMIN,YMAX))
	  hdummy[itoy]->Fill(hit_dqdx[itoy][ihit]);
      }
    }
  }
  
  //fit toy histograms and get mean value
  TF1* f = new TF1("f",Langaus,40,80,4);
  f->SetParameters(1,60,5000,5);
  f->SetParLimits(0,1,5);
  f->SetParLimits(1,50,70);
  f->SetParLimits(3,1,5);
  for(int itoy = 0; itoy < NTOYS; itoy++){
    //hdummy[itoy]->Draw();
    //gPad->Update();
    hdummy[itoy]->Fit("f","QNOLR");
    hsyst->Fill(f->GetParameter(1));
    hdummy[itoy]->Reset();
  }
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  hsyst->Reset();

  std::cout << "Global MPV = " << global_mpv << " +/- " << global_syst_error << std::endl; 
  
  //get local mean and errors
  const int nbinsz = 5;//(ZMAX-ZMIN)/STEP;
  const int nbinsy = 7;//(YMAX-YMIN)/STEP;
  double local_mpv[nbinsz][nbinsy]        = {0};
  double local_syst_error[nbinsz][nbinsy] = {0};
  double zmin_vector[nbinsz] = {0,55,110,165,220};
  double ymin_vector[nbinsy] = {0,95,195,295,395,495,595};
  std::cout << "local mean calculation" << std::endl;
  TH2F* h2 = new TH2F("h2","h2",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX);
  
  //loop over voxels
  for(int iz = 0; iz < nbinsz; iz++){
    double zmin = zmin_vector[iz];//ZMIN + iz*STEP;
    double zmax = zmin_vector[iz]+STEP;//ZMIN + (iz+1)*STEP;
    for(int iy = 0; iy < nbinsy; iy++){
      double ymin = ymin_vector[iy];//YMIN + iy*STEP;
      double ymax = ymin_vector[iy]+STEP;//YMIN + (iy+1)*STEP;
      std::cout << "bin " << iz << "/" << nbinsz << " " << iy << "/" << nbinsy << std::endl;
      std::cout << zmin << "< Z <" << zmax << " && " << ymin << " < Y <" << ymax << std::endl;
      //loop over entries
      for(int ientry = 0; ientry < nentries; ientry++){
	tree->GetEntry(ientry);
	if(ientry%20000==0)std::cout << ientry << "/" << nentries << std::endl;
	//loop over toys
	for(int itoy = 0; itoy < NTOYS; itoy++){
	  //loop over hits
	  for(int ihit = 0; ihit < NHITS; ihit++){
	    if(IsValidHit(hit_z[itoy][ihit],hit_y[itoy][ihit],zmin,zmax,ymin,ymax))
	      hdummy[itoy]->Fill(hit_dqdx[itoy][ihit]);
	  }//loop over hits
	}//loop over toys
      }//loop over entries
      //fit toy histograms and get mean value
      for(int itoy = 0; itoy < NTOYS; itoy++){
	//hdummy[itoy]->Draw();
	f->SetParameters(1,60,5000,5);
	hdummy[itoy]->Fit("f","QNOLR");
	//gPad->Update();
	hsyst->Fill(f->GetParameter(1));
	hdummy[itoy]->Reset();
      }
      //hsyst->Draw();gPad->Update();
      local_mpv[iz][iy] = hsyst->GetMean();
      local_syst_error[iz][iy] = hsyst->GetRMS();
      hsyst->Reset();
      std::cout << "Local MPV = " << local_mpv[iz][iy] << " +/- " << local_syst_error[iz][iy] << std::endl;
      h2->SetBinContent(iz+1,iy+1,sqrt(pow(global_syst_error/global_mpv,2)+pow(local_syst_error[iz][iy]/local_mpv[iz][iy],2))*100);
    }
  }//loop over voxels
  h2->Draw("colz");
  TFile* wfile = TFile::Open("yz_errormap_lifetime.root","NEW");
  h2->Write();
  wfile->Close();
}
