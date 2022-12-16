const int NMAXTRACKS = 30;
const int NMAXHITS   = 3000;
const double YMAX    = 600;
const double ZMAX    = 695;
const double YBIN    = 5;
const double ZBIN    = 5;


bool IsXContained(const double x){return (x > -360 && x < 0);}

bool IsYContained(const double y){return (y > 0 && y < YMAX);}

bool IsZContained(const double z){return (z > 0 && z < ZMAX);}

void yz(){

  //file name
  std::string filename = "/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00/localProducts_larsoft_v09_56_00_e20_prof/highland/highlandPD/dedx.root";
  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* tree  = (TTree*)rfile->Get("ana");

  //set address for variables
  int ntracks;
  float track_hit_x[NMAXTRACKS][NMAXHITS];
  float track_hit_y[NMAXTRACKS][NMAXHITS];
  float track_hit_z[NMAXTRACKS][NMAXHITS];
  float track_hit_dqdx[NMAXTRACKS][NMAXHITS];
  tree->SetBranchAddress("ntracks"       , &ntracks      );
  tree->SetBranchAddress("track_hit_x"   , track_hit_x   );
  tree->SetBranchAddress("track_hit_y"   , track_hit_y   );
  tree->SetBranchAddress("track_hit_z"   , track_hit_z   );
  tree->SetBranchAddress("track_hit_dqdx", track_hit_dqdx);

  //variables for calculations
  const int NBINSY = (int)YMAX/YBIN;
  const int NBINSZ = (int)ZMAX/ZBIN;
  std::vector<double> dqdx_values[NBINSZ][NBINSY];
  std::vector<double> all_dqdx_values;
  for(int iz = 0; iz < NBINSZ; iz++)
    for(int iy = 0; iy < NBINSY; iy++)
      dqdx_values[iz][iy].clear();
  all_dqdx_values.clear();

  //get number of entries
  int nentries = tree->GetEntries();

  TH1F* hh = new TH1F("hh","hh",9,1,10);
  int hitcounter = 0;
  int trackcounter = 0;

  //loop over entries
  for(int ientry = 0; ientry < nentries; ientry++){
    tree->GetEntry(ientry);
    if(ientry%10000 == 0)std::cout << (double)ientry/nentries*100 << " %" << std::endl;
    //loop over tracks
    for(int itrk = 0; itrk < ntracks; itrk++){
      //loop over hits
      for(int ihit = 0; ihit < NMAXHITS; ihit++){
  	//hit confinment
  	if(IsXContained(track_hit_x[itrk][ihit]) && 
  	   IsYContained(track_hit_y[itrk][ihit]) && 
  	   IsZContained(track_hit_z[itrk][ihit])){
  	  //get y z bins
  	  int z = int(track_hit_z[itrk][ihit])/ZBIN; 
  	  int y = int(track_hit_y[itrk][ihit])/YBIN;
  	  //fill vector of values
  	  dqdx_values[z][y].push_back(track_hit_dqdx[itrk][ihit]);
	  if(abs(track_hit_y[itrk][ihit]-300)<20 && abs(track_hit_z[itrk][ihit]-345)<20)
	    hitcounter++;
	}
      }
      if(hitcounter!=0)trackcounter++;
      hitcounter=0;
    }
    hh->Fill(trackcounter);
    trackcounter=0;
  }
  hh->Draw();gPad->SetLogy();gPad->Update();gPad->WaitPrimitive();std::exit(1);

  //histogram with results
  TH2F* dqdx_hist = new TH2F("dqdx_hist","dqdx_hist",NBINSZ,0,ZMAX,NBINSY,0,YMAX);
  TH2F* corr_hist = new TH2F("corr_hist","corr_hist",NBINSZ,0,ZMAX,NBINSY,0,YMAX);
  TH2F* ent_hist = new TH2F("ent_hist","ent_hist",NBINSZ,0,ZMAX,NBINSY,0,YMAX);
  TH2F* rms_hist = new TH2F("rms_hist","rms_hist",NBINSZ,0,ZMAX,NBINSY,0,YMAX);

  //loop over matrix elements
  for(int iz = 0; iz < NBINSZ; iz++){
    for(int iy = 0; iy < NBINSY; iy++){
      //loop over values per cell
      if(dqdx_values[iz][iy].size()>5){
  	for(int ival = 0; ival < (int)dqdx_values[iz][iy].size(); ival++)
  	  all_dqdx_values.push_back(dqdx_values[iz][iy][ival]);
  	//compute average dqdx per cell
  	dqdx_hist->SetBinContent(iz+1,iy+1,
				 TMath::Median(dqdx_values[iz][iy].size(),
					       &dqdx_values[iz][iy][0]));
	ent_hist->SetBinContent(iz+1,iy+1,
				dqdx_values[iz][iy].size());
	// rms_hist->SetBinContent(iz+1,iy+1,
	// 			TMath::RMS(dqdx_values[iz][iy].size(),
	// 				   &dqdx_values[iz][iy][0]) /
	// 			TMath::Median(dqdx_values[iz][iy].size(),
	// 				      &dqdx_values[iz][iy][0]) * 100);
      }
    }
  }
  
  //get global average and corrections
  double dqdx_global = TMath::Median(all_dqdx_values.size(),&all_dqdx_values[0]);
  std::cout << dqdx_global << std::endl;
  for(int iz = 0; iz < NBINSZ; iz++)
    for(int iy = 0; iy < NBINSY; iy++)
      corr_hist->SetBinContent(iz+1,iy+1, dqdx_global / dqdx_hist->GetBinContent(iz+1,iy+1));

  // rms_hist->GetZaxis()->SetRangeUser(0,100);
  // rms_hist->Draw("colz");


  //check some bins
  TH1F* h = new TH1F("h","h",100,0,100);
  TF1* f;
  // for(int iz = 0; iz < NBINSZ; iz++){
  //   for(int iy = 0; iy < NBINSY; iy++){
  //     h->Reset();
  //     //loop over values per cell
  //     if(dqdx_values[iz][iy].size()>5){
  // 	for(int i = 0; i < dqdx_values[iz][iy].size(); i++)
  // 	  h->Fill(dqdx_values[iz][iy][i]);
  //     }
  //     h->Fit("landau","Q");//Draw();gPad->Update();
  //     f = (TF1*)h->GetFunction("landau");
  //     rms_hist->SetBinContent(iz,iy,f->GetParameter(2)/f->GetParameter(1)*100);
  //   }
  // }

  // std::cout << TMath::Median(dqdx_values[50][50].size(),&dqdx_values[50][50][0]) << std::endl;
  //loop over values per cell
  if(dqdx_values[50][50].size()>5){
    for(int i = 0; i < dqdx_values[50][50].size(); i++)
      h->Fill(dqdx_values[50][50][i]);
  }
  h->Draw();
  h->Fit("landau");
  
  //rms_hist->GetZaxis()->SetRangeUser(0,100);
  //rms_hist->Draw("colz");
}
