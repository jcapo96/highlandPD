const int ZMAX = 695;
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
void toy_yz_branch_histos(){
//********************************************
  
  //file name
  std::string filename = "/dune/data/users/miagarc/toy_sce_lifetime.root";
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
  TH2F* hdummy = new TH2F("dummy","dummy",NTOYS,0,NTOYS,100,0,100);
  TH1F* hsyst = new TH1F("hsyst","hsyst",1000,0,100);

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
	  hdummy->Fill(itoy+1,hit_dqdx[itoy][ihit]);
      }
    }
  }
  
  //fit toy histograms and get mean value
  TF1* f = new TF1("f",Langaus,40,80,4);
  f->SetParameters(1,60,5000,5);
  f->SetParLimits(0,1,5);
  f->SetParLimits(1,50,70);
  f->SetParLimits(3,1,5);
  TH1D* hproj = new TH1D("hproj","hproj",100,0,100);
  for(int itoy = 0; itoy < NTOYS; itoy++){
    hproj=hdummy->ProjectionY("hproj",itoy+1,itoy+2);
    double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
    f->SetParameter(1,hmax);
    f->SetParLimits(1,hmax-5,hmax+10);
    f->SetRange(hmax-7.5,hmax+10);
    hproj->Fit("f","QNOLR");
    hsyst->Fill(f->GetParameter(1));
    hproj->Reset();
  }
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  hsyst->Reset();
  delete hdummy;

  std::cout << "Global MPV = " << global_mpv << " +/- " << global_syst_error << std::endl; 
  
  //get local mean and errors
  const int nbinsz = (ZMAX-ZMIN)/STEP;
  const int nbinsy = (YMAX-YMIN)/STEP;
  double local_mpv[nbinsz][nbinsy]        = {0};
  double local_syst_error[nbinsz][nbinsy] = {0};
  std::cout << "local mean calculation" << std::endl;
  TH2F* hdummyv[nbinsz][nbinsy];
  for(int iz = 0; iz < nbinsz; iz++){
    std::stringstream ssz;
    ssz << iz;
    for(int iy = 0; iy < nbinsy; iy++){
      std::stringstream ssy;
      ssy << iy;
      hdummyv[iz][iy] = new TH2F(("dummy_"+ssz.str()+"_"+ssy.str()+"").c_str(),
				 ("dummy_"+ssz.str()+"_"+ssy.str()+"").c_str(),
				 NTOYS,0,NTOYS,100,0,100);
    }
  }
  
  //loop over entries
  int z, y;
  for(int ientry = 0; ientry < nentries; ientry++){
    tree->GetEntry(ientry);
    if(ientry%20000==0)std::cout << ientry << "/" << nentries << std::endl;
    //loop over toys
    for(int itoy = 0; itoy < NTOYS; itoy++){
      //loop over hits
      for(int ihit = 0; ihit < NHITS; ihit++){
	if(IsValidHit(hit_z[itoy][ihit],hit_y[itoy][ihit],ZMIN,ZMAX,YMIN,YMAX)){
	  z = hit_z[itoy][ihit]/STEP;
	  y = hit_y[itoy][ihit]/STEP;
	  hdummyv[z][y]->Fill(itoy+1,hit_dqdx[itoy][ihit]);
	}
      }//loop over hits
    }//loop over toys
  }//loop over entries
    
  //fit toy histograms and get mean value
  TH2F* hf = new TH2F("hf","hf",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX);
  //loop over voxels
  for(int iz = 0; iz < nbinsz; iz++){
    for(int iy = 0; iy < nbinsy; iy++){
      //loop over toys
      for(int itoy = 0; itoy < NTOYS; itoy++){
	hproj=hdummyv[iz][iy]->ProjectionY("hproj",itoy+1,itoy+2);
	f->SetParameters(1,60,5000,5);
	double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
	f->SetParameter(1,hmax);
	f->SetParLimits(1,hmax-5,hmax+10);
	f->SetRange(hmax-7.5,hmax+10);
	hproj->Fit("f","QNOLR");
	hsyst->Fill(f->GetParameter(1));
	hproj->Reset();
      }
      local_mpv[iz][iy] = hsyst->GetMean();
      local_syst_error[iz][iy] = hsyst->GetRMS();
      hsyst->Reset();
      delete hdummyv[iz][iy];
      std::cout << iz << " " << iy << " Local MPV = " << local_mpv[iz][iy] << " +/- " << local_syst_error[iz][iy] << std::endl;
      hf->SetBinContent(iz+1,iy+1,sqrt(pow(global_syst_error/global_mpv,2)+pow(local_syst_error[iz][iy]/local_mpv[iz][iy],2))*100);
    }
  }
    
  hf->Draw("colz");
  TFile* wfile = TFile::Open("yz_errormap_complete.root","NEW");
  hf->Write();
  wfile->Close();
}
