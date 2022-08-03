#define Calorimetry_C

#include "Calorimetry.hxx"
#include "Parameters.hxx"
#include "SpaceCharge.hxx"

//********************************************************************
Calorimetry::Calorimetry(){
//********************************************************************

  _h_norm = NULL;
  for(int i = 0; i < 3; i++)_h_Xcorr[i] = NULL;
  for(int i = 0; i < 3; i++)for(int j = 0; j < 2; j++)_h_YZcorr[i][j] = NULL;
  
  _Efield  = -999.;
  _ModBoxA = -999.;
  _ModBoxB = -999.;
  for(int i = 0; i < 3; i++)_CalAreaConstants[i] = -999.;
}

//********************************************************************
Calorimetry::~Calorimetry(){
//********************************************************************

  delete _h_norm;
  for(int i = 0; i < 3; i++)delete _h_Xcorr[i];
  for(int i = 0; i < 3; i++)for(int j = 0; j < 2; j++)delete _h_YZcorr[i][j];

}

//********************************************************************
void Calorimetry::Initialize(){
//********************************************************************

  //read data files and create correction histograms
  CreateNormHistogram();
  CreateXcorrHistogram();
  CreateYZcorrHistogram();

  _Efield  = 0.4867; //probably all these values should be added in the parameters file
  _ModBoxA = 0.930;
  _ModBoxB = 0.212; 

  _CalAreaConstants[0] = 1.04e-3;
  _CalAreaConstants[1] = 1.04e-3;
  _CalAreaConstants[2] = 1.0156e-3;

}

//********************************************************************
void Calorimetry::CreateNormHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dEdx_norm.txt";
  FILE* tf = fopen(filename.c_str(), "r");
  
  //check it exits
  if(!tf){
    std::cout << "Can't open file " << filename << std::endl;
    exit(1);
  }

  //normalization correction has a value for each plane
  _h_norm = new TH1F("hNorm","hNorm",3,0,3);

  //ignore header
  char line[100];
  fgets(line,100,tf);
  
  //get values
  double bin,norm,error,dummy;
  while(fscanf(tf, "%lf,%lf,%lf,%lf", &bin, &dummy, &norm, &error) == 4){
    _h_norm->SetBinContent(bin,norm);
    _h_norm->SetBinError(bin,error);
  }

  fclose(tf);
}

//********************************************************************
void Calorimetry::CreateXcorrHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dEdx_Xcorr.txt";
  FILE* tf = fopen(filename.c_str(), "r");
  
  //check it exits
  if(!tf){
    std::cout << "Can't open file " << filename << std::endl;
    exit(1);
  }

  //ignore header
  char line[100];
  fgets(line,100,tf);
  
  //get values
  std::vector<double> x_vector[3], corr_vector[3], error_vector[3];
  for(int i = 0; i < 3; i++){
    x_vector[i].clear();
    corr_vector[i].clear();
    error_vector[i].clear();
  }
  double bin, corr, error, x, dummy1, dummy2;
  int plane;
  while(fscanf(tf, "%lf,%lf,%lf,%lf,%lf,%lf", &bin, &dummy1, &corr, &error, &x, &dummy2) == 6){
    plane = (int)bin/10000;
    x_vector[plane].push_back(x);    
    corr_vector[plane].push_back(corr);    
    error_vector[plane].push_back(error);    
  }

  //create histograms
  int nbins    = x_vector[0].size();
  double xmax  = x_vector[0].back();
  double xmin  = x_vector[0][0];
  double width = (xmax-xmin)/(nbins-1);

  for(int iplane = 0; iplane < 3; iplane++){
    std::stringstream ssi;
    ssi << iplane;
    _h_Xcorr[iplane] = new TH1F(("hXcorr_plane"+ssi.str()+"").c_str(),
				("hXcorr_plane"+ssi.str()+"").c_str(),
				nbins,xmin-width/2,xmax+width/2);
    for(int ibin = 0; ibin < (int)x_vector[iplane].size(); ibin++){
      _h_Xcorr[iplane]->SetBinContent(ibin+1,corr_vector[iplane][ibin]);
      _h_Xcorr[iplane]->SetBinError(ibin+1,error_vector[iplane][ibin]);
    }
  }    
  fclose(tf);
}

//********************************************************************
void Calorimetry::CreateYZcorrHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dEdx_YZcorr.txt";
  FILE* tf = fopen(filename.c_str(), "r");
  
  //check it exits
  if(!tf){
    std::cout << "Can't open file " << filename << std::endl;
    exit(1);
  }

  //ignore header
  char line[100];
  fgets(line,100,tf);

  //get histogram limits
  double ymin, ymax, zmin, zmax, ywidth, zwidth;
  ymin = zmin = ywidth = zwidth = 100;
  ymax = zmax = -100;

  double channel, corr, error, y, z, dummy;
  while(fscanf(tf, "%lf,%lf,%lf,%lf,%lf,%lf", &channel, &dummy, &corr, &error, &y, &z) == 6){
    if(y < ymin)ymin = y;
    if(y > ymax)ymax = y;
    if(z < zmin)zmin = z;
    if(z > zmax)zmax = z;
    if(ymax - ymin != 0 && ymax - ymin < ywidth)ywidth = ymax - ymin;
    if(zmax - zmin != 0 && zmax - zmin < zwidth)zwidth = zmax - zmin;
  }
  fclose(tf);

  //create histograms
  int nbinsy, nbinsz;
  nbinsy = (ymax-ymin)/ywidth + 1;
  nbinsz = (zmax-zmin)/zwidth + 1;
  
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      std::stringstream ssi, ssj;
      ssi << i;
      ssj << j;
      _h_YZcorr[i][j] = new TH2F(("hYZcorr_plane"+ssi.str()+"_side"+ssj.str()+"").c_str(),
				 ("hYZcorr_plane"+ssi.str()+"_side"+ssj.str()+"").c_str(),
				 nbinsz, zmin-zwidth/2, zmax+zwidth/2,
				 nbinsy, ymin-ywidth/2, ymax+ywidth/2);
    }
  }

  tf = fopen(filename.c_str(), "r");
  fgets(line,100,tf);
  int plane, side;
  while(fscanf(tf, "%lf,%lf,%lf,%lf,%lf,%lf", &channel, &dummy, &corr, &error, &y, &z) == 6){
    plane = channel/10000000;
    side = (channel - plane*10000000)/1000000;
    _h_YZcorr[plane][side]->Fill(z,y,corr);
  }
  fclose(tf);    
}


//********************************************************************
void Calorimetry::CalibrateHit(AnaHitPD &hit) const {
//********************************************************************
  
  //check plane ID
  if(hit.PlaneID < 0 || hit.PlaneID > 2){
    std::cout << "wrong plane!" << std::endl;
    return;
  }

  //initialize dEdx value to dQdx value
  hit.dQdx = hit.dQdx_NoSCE;

    //normalization correction
  //ApplyNormCorrection(hit);

  //X correction
  ApplyXCorrection(hit);

  //YZ correction
  ApplyYZCorrection(hit);

  //recombination
  ApplyRecombination(hit);  
}

//********************************************************************
void Calorimetry::ApplyNormCorrection(AnaHitPD &hit) const {
//********************************************************************

  if(!_h_norm){
    std::cout << "normalization histogram does not exist, you should call Initialize method!" << std::endl;
    return;
  }

  hit.dQdx *= _h_norm->GetBinContent(hit.PlaneID); 

}

//********************************************************************
void Calorimetry::ApplyXCorrection(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    if(!_h_Xcorr[i]){
      std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
      return;
    } 
  }

  int bin = _h_Xcorr[hit.PlaneID]->FindBin(hit.Position.X());
  hit.dQdx *= _h_Xcorr[hit.PlaneID]->GetBinContent(bin);
}

//********************************************************************
void Calorimetry::ApplyYZCorrection(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      if(!_h_YZcorr[i][j]){
	std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
	return;
      } 
    }
  }
  
  int side = (int)(hit.Position.X() > 0); 
  int bin  = _h_YZcorr[hit.PlaneID][side]->FindBin(hit.Position.Z(),hit.Position.Y());
  hit.dQdx *= _h_YZcorr[hit.PlaneID][side]->GetBinContent(bin);
}

//********************************************************************
void Calorimetry::ApplyRecombination(AnaHitPD &hit) const {
//********************************************************************

  double dQdx_e = GetElectronsFromADCArea(hit.dQdx,hit.PlaneID); //electrons from ADC area
  double rho = 1.38877;                                          //LAr density in g/cm³
  double Wion = 1000./4.237e7;                                   //Wion in MeV/e
  double Efield = _Efield;                                       //kV/cm

  //E field SCE
  SpaceCharge* sce = new SpaceCharge();
  sce->Initialize();
  //the -1 is added ad hoc. I don't know the reason but there is a -1 difference wrt the larsoft calculation
  TVector3 offset = -1*(sce->GetCalEfieldOffsets(hit.Position,hit.TPCid)); 
  Efield = sqrt(pow(Efield*(1+offset.X()),2) + pow(Efield*offset.Y(),2) + pow(Efield*offset.Z(),2));

  //calculate recombination factors
  double Beta = _ModBoxB / (rho * Efield);
  double Alpha = _ModBoxA;
  hit.dEdx = (exp(Beta * Wion * dQdx_e) - Alpha) / Beta;
}
