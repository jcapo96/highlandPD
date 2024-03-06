//*********************************************************************************
double GetWeight(TH3F*h1, int x1, int y1, int z1, TH3F*h2, int x2, int y2, int z2){
//*********************************************************************************

  //get voxel borders
  double xmin1 = h1->GetXaxis()->GetBinCenter(x1)-h1->GetXaxis()->GetBinWidth(x1)/2;
  double xmax1 = h1->GetXaxis()->GetBinCenter(x1)+h1->GetXaxis()->GetBinWidth(x1)/2;
  double ymin1 = h1->GetYaxis()->GetBinCenter(y1)-h1->GetYaxis()->GetBinWidth(y1)/2;
  double ymax1 = h1->GetYaxis()->GetBinCenter(y1)+h1->GetYaxis()->GetBinWidth(y1)/2;
  double zmin1 = h1->GetZaxis()->GetBinCenter(z1)-h1->GetZaxis()->GetBinWidth(z1)/2;
  double zmax1 = h1->GetZaxis()->GetBinCenter(z1)+h1->GetZaxis()->GetBinWidth(z1)/2;

  //get voxel borders
  double xmin2 = h2->GetXaxis()->GetBinCenter(x2)-h2->GetXaxis()->GetBinWidth(x2)/2;
  double xmax2 = h2->GetXaxis()->GetBinCenter(x2)+h2->GetXaxis()->GetBinWidth(x2)/2;
  double ymin2 = h2->GetYaxis()->GetBinCenter(y2)-h2->GetYaxis()->GetBinWidth(y2)/2;
  double ymax2 = h2->GetYaxis()->GetBinCenter(y2)+h2->GetYaxis()->GetBinWidth(y2)/2;
  double zmin2 = h2->GetZaxis()->GetBinCenter(z2)-h2->GetZaxis()->GetBinWidth(z2)/2;
  double zmax2 = h2->GetZaxis()->GetBinCenter(z2)+h2->GetZaxis()->GetBinWidth(z2)/2;

  double xa,xb,ya,yb,za,zb;
  if(xmin1>xmin2)xa=xmin1;
  else           xa=xmin2;
  if(xmax1<xmax2)xb=xmax1;
  else           xb=xmax2;
  if(ymin1>ymin2)ya=ymin1;
  else           ya=ymin2;
  if(ymax1<ymax2)yb=ymax1;
  else           yb=ymax2;
  if(zmin1>zmin2)za=zmin1;
  else           za=zmin2;
  if(zmax1<zmax2)zb=zmax1;
  else           zb=zmax2;

  //std::cout << (xb-xa)*(yb-ya)*(zb-za) / ((xmax2-xmin2)*(ymax2-ymin2)*(zmax2-zmin2)) << std::endl;
  return (xb-xa)*(yb-ya)*(zb-za) / ((xmax2-xmin2)*(ymax2-ymin2)*(zmax2-zmin2));

}

//*********************************************************************************
double GetAverageValue(TH3F*h1, TH3F*h2, int x2, int y2, int z2){
//*********************************************************************************

  //get voxel borders
  double xmin2 = h2->GetXaxis()->GetBinCenter(x2)-h2->GetXaxis()->GetBinWidth(x2)/2;
  double xmax2 = h2->GetXaxis()->GetBinCenter(x2)+h2->GetXaxis()->GetBinWidth(x2)/2;
  double ymin2 = h2->GetYaxis()->GetBinCenter(y2)-h2->GetYaxis()->GetBinWidth(y2)/2;
  double ymax2 = h2->GetYaxis()->GetBinCenter(y2)+h2->GetYaxis()->GetBinWidth(y2)/2;
  double zmin2 = h2->GetZaxis()->GetBinCenter(z2)-h2->GetZaxis()->GetBinWidth(z2)/2;
  double zmax2 = h2->GetZaxis()->GetBinCenter(z2)+h2->GetZaxis()->GetBinWidth(z2)/2;

  double average = 0;
  int nvoxels = 0;
  double totalweight = 0;

  for(int x = 1; x <= h1->GetNbinsX(); x++){
    double xmin = h1->GetXaxis()->GetBinCenter(x)-h1->GetXaxis()->GetBinWidth(x)/2;
    double xmax = h1->GetXaxis()->GetBinCenter(x)+h1->GetXaxis()->GetBinWidth(x)/2;
    if((xmin < xmax2 && xmin > xmin2) || (xmax < xmax2 && xmax > xmin2)){
      for(int y = 1; y <= h1->GetNbinsY(); y++){
	double ymin = h1->GetYaxis()->GetBinCenter(y)-h1->GetYaxis()->GetBinWidth(y)/2;
	double ymax = h1->GetYaxis()->GetBinCenter(y)+h1->GetYaxis()->GetBinWidth(y)/2;
	if((ymin < ymax2 && ymin > ymin2) || (ymax < ymax2 && ymax > ymin2)){
	  for(int z = 1; z <= h1->GetNbinsZ(); z++){
	    double zmin = h1->GetZaxis()->GetBinCenter(z)-h1->GetZaxis()->GetBinWidth(z)/2;
	    double zmax = h1->GetZaxis()->GetBinCenter(z)+h1->GetZaxis()->GetBinWidth(z)/2;
	    if((zmin < zmax2 && zmin > zmin2) || (zmax < zmax2 && zmax > zmin2)){
	      double weight = GetWeight(h1,x,y,z,h2,x2,y2,z2);
	      totalweight += weight;
	      average += h1->GetBinContent(x,y,z)*weight;
	      nvoxels++;
	    }
	  }
	}
      }
    }
  }
  if(totalweight==1)return average;
  else return -99999;
  //return average;///nvoxels;
}

//*********************************************************************************
std::vector<int> GetClosestVoxel(TH3F* h1, TH3F* h2, int x2, int y2, int z2){
//*********************************************************************************

  std::vector<int> bins(3);
  double dif = 999;
  TVector3 center1(0,0,0);
  TVector3 center2(h2->GetXaxis()->GetBinCenter(x2),
		   h2->GetYaxis()->GetBinCenter(y2),
		   h2->GetZaxis()->GetBinCenter(z2));

  for(int x = 1; x <= h1->GetNbinsX(); x++){
    for(int y = 1; y <= h1->GetNbinsY(); y++){
      for(int z = 1; z <= h1->GetNbinsZ(); z++){
	center1.SetXYZ(h1->GetXaxis()->GetBinCenter(x),
		       h1->GetYaxis()->GetBinCenter(y),
		       h1->GetZaxis()->GetBinCenter(z));
	TVector3 distance = center2-center1;
	if(distance.Mag()<dif){
	  dif = distance.Mag();
	  bins[0] = x;
	  bins[1] = y;
	  bins[2] = z;
	}
      }
    }
  }
  center1.SetXYZ(h1->GetXaxis()->GetBinCenter(bins[0]),
		 h1->GetYaxis()->GetBinCenter(bins[1]),
		 h1->GetZaxis()->GetBinCenter(bins[2]));
  
  center2.Print();
  center1.Print();
  return bins;
}

//*********************************************************************************
void maps(){
//*********************************************************************************

  //nominal sce map
  std::string filename1 = "/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00/localProducts_larsoft_v09_56_00_e20_prof/highland/highlandPD/src/pdUtils/data/SCE_DataDriven_180kV_v4.root";
  //open root file
  TFile* rfile1 = TFile::Open(filename1.c_str());

  //alternative sce map
  std::string filename2 = "/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00/localProducts_larsoft_v09_56_00_e20_prof/highland/highlandPD/src/pdUtils/data/SCE_Alternate_v4_EField.root";
  //open root file
  TFile* rfile2 = TFile::Open(filename2.c_str());

  //get maps
  TH3F* h1 = (TH3F*)rfile1->Get("RecoBkwd_Displacement_Y_Neg");
  TH3F* h2 = (TH3F*)rfile2->Get("RecoBkwd_Displacement_Y_Neg");

  //histogram for difference
  TH1F* h = new TH1F("h","h",20,-10,10);

  for(int x = 1; x <= h2->GetNbinsX(); x++){
    for(int y = 1; y <= h2->GetNbinsY(); y++){
      for(int z = 1; z <= h2->GetNbinsZ(); z++){
  	//look for closest original bin
  	//std::vector<int> bins = GetClosestVoxel(h1,h2,x,y,z);
  	float altval = h2->GetBinContent(x,y,z);
  	float orival = GetAverageValue(h1,h2,x,y,z);//h1->GetBinContent(bins[0],bins[1],bins[2]);
	if(altval==0 || orival==-99999)continue;
	std::cout << altval << " " << orival << std::endl;
  	h->Fill((altval-orival));
  	//if(abs((altval-orival) / orival *100)>10)
  	//std::cout << "bin (" << x << "," << y << "," << z << ")" << std::endl;
      }
    }
  }
   
  h->Draw();

  // std::cout << h1->GetXaxis()->GetBinWidth(4) << std::endl;
  // std::cout << h1->GetYaxis()->GetBinWidth(4) << std::endl;
  // std::cout << h1->GetZaxis()->GetBinWidth(4) << std::endl;
}
