void binning(){

  std::string filename = "/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00/localProducts_larsoft_v09_56_00_e20_prof/highland/highlandPD/src/pdUtils/data/SCE_DataDriven_180kV_v4.root";
  TFile* rfile = TFile::Open(filename.c_str());
  
  TH3F* h = (TH3F*)rfile->Get("RecoBkwd_Displacement_X_Neg");

  int nbinsx = h->GetNbinsX();
  int nbinsy = h->GetNbinsY();
  int nbinsz = h->GetNbinsZ();

  double xmin,xmax,ymin,ymax,zmin,zmax;
  for(int x = 0; x < nbinsx; x++){
    xmin = h->GetXaxis()->GetBinCenter(x+1) - 0.5*h->GetXaxis()->GetBinWidth(x+1);
    xmax = h->GetXaxis()->GetBinCenter(x+1) + 0.5*h->GetXaxis()->GetBinWidth(x+1);
    for(int y = 0; y < nbinsy; y++){
      ymin = h->GetYaxis()->GetBinCenter(y+1) - 0.5*h->GetYaxis()->GetBinWidth(y+1);
      ymax = h->GetYaxis()->GetBinCenter(y+1) + 0.5*h->GetYaxis()->GetBinWidth(y+1);
      for(int z = 0; z < nbinsz; z++){
  	zmin = h->GetZaxis()->GetBinCenter(z+1) - 0.5*h->GetZaxis()->GetBinWidth(z+1);
  	zmax = h->GetZaxis()->GetBinCenter(z+1) + 0.5*h->GetZaxis()->GetBinWidth(z+1);
  	std::cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << " " << 0.05 << std::endl;
      }
    }
  }
}
