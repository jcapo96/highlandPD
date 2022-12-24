const int XMAX = 0;
const int XMIN = -360;
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

TFile* ef = TFile::Open("/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00/localProducts_larsoft_v09_56_00_e20_prof/highland/highlandPD/src/pdUtils/data/SCE_DataDriven_180kV_v4.root");
TH3F* Ex  =(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F* Ey  =(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F* Ez  =(TH3F*)ef->Get("Reco_ElecField_Z_Neg");

float ElectricField(float x,float y,float z){
  float E0 = 0.4867;
  float ex = E0+E0*Ex->GetBinContent(Ex->FindBin(x,y,z));
  float ey = 0.0+E0*Ey->GetBinContent(Ex->FindBin(x,y,z));
  float ez = 0.0+E0*Ex->GetBinContent(Ex->FindBin(x,y,z));
  return sqrt(ex*ex+ey*ey+ez*ez);
}

//********************************************
double Recombination(const double E){
//********************************************

  float LAr_density=1.39;
  float alp=0.93;
  float bet=0.212;
  float dedx=1.9;
  float xsi=bet*dedx/(LAr_density*E);
  float xsi0=bet*dedx/(LAr_density*0.4867);
  float rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}

//********************************************
bool IsValidHit(const double x, const double y, const double z,
		const double xmin, const double xmax){
//********************************************

  return (x > xmin && x < xmax && 
	  y > 0    && y < 600  && 
	  z > 0    && z < 250);
}

//********************************************
void toy_x_branch(){
//********************************************
  
  //file name
  std::string filename = "/dune/data/users/miagarc/toy.root";
  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* tree  = (TTree*)rfile->Get("all_syst");

  //tree variables
  Float_t hit_x [NTOYS][NHITS]    = {0};
  Float_t hit_y [NTOYS][NHITS]    = {0};
  Float_t hit_z [NTOYS][NHITS]    = {0};
  Float_t hit_dqdx [NTOYS][NHITS] = {0};

  tree->SetBranchAddress("toy_hit_x"   , hit_x   );
  tree->SetBranchAddress("toy_hit_y"   , hit_y   );
  tree->SetBranchAddress("toy_hit_z"   , hit_z   );
  tree->SetBranchAddress("toy_hit_dqdx", hit_dqdx);

  int nentries = tree->GetEntries();
  
  //get yz correction map
  TFile* yzfile = TFile::Open("YZcorrection.root");
  TH2F* yzhisto = (TH2F*)yzfile->Get("corr_hist");

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
    if(ientry%10000==0)std::cout << ientry << "/" << nentries << std::endl;
    //loop over toys
    for(int itoy = 0; itoy < NTOYS; itoy++){
      //loop over hits
      for(int ihit = 0; ihit < NHITS; ihit++){
	if(IsValidHit(hit_x[itoy][ihit],hit_y[itoy][ihit],hit_z[itoy][ihit],XMIN,XMAX)){
	  //apply corrections
	  double yzcorr = yzhisto->GetBinContent(yzhisto->GetXaxis()->FindBin(hit_z[itoy][ihit]),
						 yzhisto->GetYaxis()->FindBin(hit_y[itoy][ihit]));
	  double refact = Recombination(ElectricField(hit_x[itoy][ihit],
						      hit_y[itoy][ihit],
						      hit_z[itoy][ihit]));
	  hdummy[itoy]->Fill(hit_dqdx[itoy][ihit]*yzcorr*refact);
	}
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
    std::cout << "toy ploy " << itoy << std::endl;
    //hdummy[itoy]->Draw();
    hdummy[itoy]->Fit("f","QNOLR");
    hsyst->Fill(f->GetParameter(1));
    hdummy[itoy]->Reset();
  }
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  hsyst->Reset();

  std::cout << "Global MPV = " << global_mpv << " +/- " << global_syst_error << std::endl; 
  
  //get local mean and errors
  const int nbinsx = 11;//(XMAX-XMIN)/STEP;
  double local_mpv[nbinsx] = {0};
  double local_syst_error[nbinsx] = {0};
  double xmin_vector[nbinsx] = {-360,-320,-280,-240,-200,-160,-120,-80,-40,-20,-5};
  std::cout << "local mean calculation" << std::endl;
  TH1F* h1f = new TH1F("h1f","h1f",nbinsx,XMIN,XMAX);
  
  //loop over voxels
  for(int ix = 0; ix < nbinsx; ix++){
    double xmin = xmin_vector[ix];//XMIN + ix*STEP;
    double xmax = xmin_vector[ix]+STEP;//XMIN + (ix+1)*STEP;
    std::cout << "bin " << ix << "/" << nbinsx << std::endl;
    //loop over entries
    for(int ientry = 0; ientry < nentries; ientry++){
      tree->GetEntry(ientry);
      if(ientry%20000==0)std::cout << ientry << "/" << nentries << std::endl;
      //loop over toys
      for(int itoy = 0; itoy < NTOYS; itoy++){
	//loop over hits
	for(int ihit = 0; ihit < NHITS; ihit++){
	  if(IsValidHit(hit_x[itoy][ihit],hit_y[itoy][ihit],hit_z[itoy][ihit],xmin,xmax)){
	    //apply corrections
	    double yzcorr = yzhisto->GetBinContent(yzhisto->GetXaxis()->FindBin(hit_z[itoy][ihit]),
						   yzhisto->GetYaxis()->FindBin(hit_y[itoy][ihit]));
	    double refact = Recombination(ElectricField(hit_x[itoy][ihit],
							hit_y[itoy][ihit],
							hit_z[itoy][ihit]));
	    hdummy[itoy]->Fill(hit_dqdx[itoy][ihit]*yzcorr*refact);
	  }
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
    local_mpv[ix] = hsyst->GetMean();
    local_syst_error[ix] = hsyst->GetRMS();
    hsyst->Reset();
    std::cout << "Local MPV = " << local_mpv[ix] << " +/- " << local_syst_error[ix] << std::endl;
    h1f->SetBinContent(ix+1,sqrt(pow(global_syst_error/global_mpv,2)+pow(local_syst_error[ix]/local_mpv[ix],2))*100);
  }//loop over voxels
  h1f->Draw();
  TFile* wfile = TFile::Open("x_errormap.root","NEW");
  h1f->Write();
  wfile->Close();
}
