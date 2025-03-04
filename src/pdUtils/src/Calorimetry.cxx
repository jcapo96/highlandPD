#define Calorimetry_C

#include "Calorimetry.hxx"
#include "Parameters.hxx"

//********************************************************************
Calorimetry::Calorimetry(){
//********************************************************************

  _h_norm = NULL;
  for(int i = 0; i < 3; i++)_h_XCal[i] = NULL;
  for(int i = 0; i < 3; i++)for(int j = 0; j < 2; j++)_h_YZCal[i][j] = NULL;
  
  _Efield  = -999.;
  _ModBoxA = -999.;
  _ModBoxB = -999.;

  _Lifetime = -999.;
  _Vdrift   = -999.;
  _APAx     = -999.;

  for(int i = 0; i < 3; i++)_CalAreaConstants[i] = -999.;

  _sce = NULL;
}

//********************************************************************
Calorimetry::~Calorimetry(){
//********************************************************************

  delete _h_norm;
  for(int i = 0; i < 3; i++)delete _h_XCal[i];
  for(int i = 0; i < 3; i++)for(int j = 0; j < 2; j++)delete _h_YZCal[i][j];

  delete _sce;
}

//********************************************************************
void Calorimetry::SetSCE(SpaceCharge* sce, bool remove){
//********************************************************************

  if(remove)
    delete _sce;
  _sce = sce;
}

//********************************************************************
void Calorimetry::Initialize(){
//********************************************************************

  //read data files and create correction histograms
  CreateNormHistogram();
  CreateXCalHistogram();
  CreateYZCalHistogram();

  _Efield  = 0.4867; //probably all these values should be added in the parameters file
  _ModBoxA = 0.930;
  _ModBoxB = 0.212; 

  _Lifetime = ND::params().GetParameterI("pdUtils.Calorimetry.ElectronLifetime.MC"); //initialize to MC by default
  _Vdrift   = 0.156461;
  _APAx     = 368.351;

  _CalAreaConstants[0] = 1.04e-3;
  _CalAreaConstants[1] = 1.04e-3;
  _CalAreaConstants[2] = 1.0156e-3;

  if(!_sce){
    _sce = new SpaceCharge();
    _sce->Initialize();
  }
}

//********************************************************************
void Calorimetry::Initialize(SpaceCharge* sce){
//********************************************************************

  _sce = sce;
  Initialize();
}

//********************************************************************
void Calorimetry::ResetModBoxParameters(){
//********************************************************************

  _ModBoxA = 0.930; //should be done by reading parameters file
  _ModBoxB = 0.212; 
}

//********************************************************************
void Calorimetry::CreateNormHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dQdx_norm.txt";
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
void Calorimetry::CreateXCalHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dQdx_XCal.txt";
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
    _h_XCal[iplane] = new TH1F(("hXCal_plane"+ssi.str()+"").c_str(),
				("hXCal_plane"+ssi.str()+"").c_str(),
				nbins,xmin-width/2,xmax+width/2);
    for(int ibin = 0; ibin < (int)x_vector[iplane].size(); ibin++){
      _h_XCal[iplane]->SetBinContent(ibin+1,corr_vector[iplane][ibin]);
      _h_XCal[iplane]->SetBinError(ibin+1,error_vector[iplane][ibin]);
    }
  }    
  fclose(tf);
}

//********************************************************************
void Calorimetry::CreateYZCalHistogram(){
//********************************************************************

  std::string filename = std::string(getenv("PDUTILSROOT"))+"/data/mc_prod4_dQdx_YZCal.txt";
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
      _h_YZCal[i][j] = new TH2F(("hYZCal_plane"+ssi.str()+"_side"+ssj.str()+"").c_str(),
				 ("hYZCal_plane"+ssi.str()+"_side"+ssj.str()+"").c_str(),
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
    _h_YZCal[plane][side]->Fill(z,y,corr);
  }
  fclose(tf);    
}

//********************************************************************
void Calorimetry::SetLifetime(const AnaEventPD &event){
//********************************************************************
  
  int run = event.EventInfo->Run;
  //if(event.GetIsMC())
  if(run>10000)
    _Lifetime = ND::params().GetParameterI("pdUtils.Calorimetry.ElectronLifetime.MC");
  else{
    std::stringstream ssrun;
    ssrun << run;
    _Lifetime = ND::params().GetParameterI(("pdUtils.Calorimetry.ElectronLifetime.Run"+ssrun.str()+"").c_str());
  }
}

//********************************************************************
void Calorimetry::CalibratedQdx(AnaParticlePD* part) const {
//********************************************************************
  
  if(!part)return;
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    CalibratedQdx(part->Hits[2][ihit]);
}

//********************************************************************
void Calorimetry::CalibratedQdx(AnaHitPD &hit) const {
//********************************************************************
  
  //check plane ID
  if(hit.PlaneID < 0 || hit.PlaneID > 2){
    std::cout << "wrong plane!" << std::endl;
    return;
  }

  //initialize dQdx value to dQdx non corrrected value value
  hit.dQdx = hit.dQdx_elife;

  //normalization correction
  //ApplyNormalization(hit);

  //X correction
  ApplyXCalibration(hit);

  //YZ correction
  ApplyYZCalibration(hit);
}

//********************************************************************
void Calorimetry::ApplyNormalization(AnaHitPD &hit) const {
//********************************************************************

  if(!_h_norm){
    std::cout << "normalization histogram does not exist, you should call Initialize method!" << std::endl;
    return;
  }

  hit.dQdx *= _h_norm->GetBinContent(hit.PlaneID); 
}

//********************************************************************
void Calorimetry::ApplyXCalibration(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    if(!_h_XCal[i]){
      std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
      return;
    } 
  }

  int bin = _h_XCal[hit.PlaneID]->FindBin(hit.Position.X());
  hit.dQdx *= _h_XCal[hit.PlaneID]->GetBinContent(bin);
}

//********************************************************************
void Calorimetry::ApplyYZCalibration(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      if(!_h_YZCal[i][j]){
	std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
	return;
      } 
    }
  }
  
  int side = (int)(hit.Position.X() > 0); 
  int bin  = _h_YZCal[hit.PlaneID][side]->FindBin(hit.Position.Z(),hit.Position.Y());
  hit.dQdx *= _h_YZCal[hit.PlaneID][side]->GetBinContent(bin);
}

//********************************************************************
void Calorimetry::ApplyRecombination(AnaHitPD &hit) const {
//********************************************************************

  double dQdx_e = GetElectronsFromADCArea(hit.dQdx,hit.PlaneID); //electrons from ADC area
  double rho = 1.38877;                                          //LAr density in g/cm³
  double Wion = 1000./4.237e7;                                   //Wion in MeV/e
  double Efield = _Efield;                                       //kV/cm

  //the -1 is added ad hoc. I don't know the reason but there is a -1 difference wrt the larsoft calculation
  TVector3 offset = -1*(_sce->GetCalEfieldOffsets(hit.Position,hit.TPCid)); 
  Efield = sqrt(pow(Efield*(1+offset.X()),2) + pow(Efield*offset.Y(),2) + pow(Efield*offset.Z(),2));

  //calculate recombination factors
  double Beta = _ModBoxB / (rho * Efield);
  double Alpha = _ModBoxA;
  hit.dEdx = (exp(Beta * Wion * dQdx_e) - Alpha) / Beta;
}

//********************************************************************
void Calorimetry::ApplyRecombination(AnaParticlePD* part) const {
//********************************************************************
  
  if(!part)return;
  if(part->Hits[2].empty())return;
  
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    ApplyRecombination(part->Hits[2][ihit]);

}

//********************************************************************
double Calorimetry::GetNormalization(AnaHitPD &hit) const {
//********************************************************************

  if(!_h_norm){
    std::cout << "normalization histogram does not exist, you should call Initialize method!" << std::endl;
    return 0.;
  }

  return _h_norm->GetBinContent(hit.PlaneID); 
}

//********************************************************************
double Calorimetry::GetXCalibration(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    if(!_h_XCal[i]){
      std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
      return 0.;
    } 
  }

  int bin = _h_XCal[hit.PlaneID]->FindBin(hit.Position.X());
  return _h_XCal[hit.PlaneID]->GetBinContent(bin);
}

//********************************************************************
double Calorimetry::GetYZCalibration(AnaHitPD &hit) const {
//********************************************************************

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      if(!_h_YZCal[i][j]){
	std::cout << "X correction histogram does not exist, you should call Initialize method!" << std::endl;
	return 0.;
      } 
    }
  }
  
  int side = (int)(hit.Position.X() > 0); 
  int bin  = _h_YZCal[hit.PlaneID][side]->FindBin(hit.Position.Z(),hit.Position.Y());
  return _h_YZCal[hit.PlaneID][side]->GetBinContent(bin);
}

//********************************************************************
void Calorimetry::ApplySCECorrection(AnaParticlePD* part) const {
//********************************************************************

  if(!part)return;
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    ApplySCECorrection(part->Hits[2][ihit]);
}


//********************************************************************
void Calorimetry::ApplySCECorrection(AnaHitPD &hit) const {
//********************************************************************

  double charge = hit.dQdx_NoSCE * hit.Pitch_NoSCE;
  TVector3 hitpos = hit.Position_NoSCE;
  TVector3 hitdir = hit.Direction_NoSCE;
  
  //compute projection of YZ plane to wire width (Z direction for collection)
  double AngleToVert = 0; //only for collection plane, which is vertical. 
  double cosgamma = abs(sin(AngleToVert)*hitdir.Y() + cos(AngleToVert)*hitdir.Z());
  double pitch = 0.4792 / cosgamma; //collection wire pitch
  
  //correct pitch by SCE effect
  TVector3 dirProjection(hitpos.X()+pitch*hitdir.X(),hitpos.Y()+pitch*hitdir.Y(),hitpos.Z()+pitch*hitdir.Z());
  TVector3 dirOffset = _sce->GetCalPosOffsets(dirProjection, hit.TPCid);
  TVector3 posOffset = _sce->GetCalPosOffsets(hitpos, hit.TPCid);
  TVector3 dirCorrection(pitch*hitdir.X() - dirOffset.X() + posOffset.X(),
			 pitch*hitdir.Y() + dirOffset.Y() - posOffset.Y(),
			 pitch*hitdir.Z() + dirOffset.Z() - posOffset.Z());
  pitch = dirCorrection.Mag();
  hit.Pitch = pitch;
  hit.dQdx_SCE = charge / pitch;
}

//********************************************************************
void Calorimetry::ApplyLifetimeCorrection(AnaParticlePD* part) const {
//********************************************************************

  if(!part)return;
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    ApplyLifetimeCorrection(part->Hits[2][ihit]);
}


//********************************************************************
void Calorimetry::ApplyLifetimeCorrection(AnaHitPD &hit) const {
//********************************************************************
  
  double xcorr = exp((_APAx-abs(hit.Position.X()))/(_Lifetime*_Vdrift));
  hit.dQdx_elife = hit.dQdx_SCE*xcorr;
}

//********************************************************************
void Calorimetry::UndoLifetimeCorrection(AnaParticlePD* part) const {
//********************************************************************

  if(!part)return;
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    UndoLifetimeCorrection(part->Hits[2][ihit]);
}


//********************************************************************
void Calorimetry::UndoLifetimeCorrection(AnaHitPD &hit) const {
//********************************************************************
  
  double xcorr = exp((_APAx-abs(hit.Position.X()))/(_Lifetime*_Vdrift));
  hit.dQdx_SCE = hit.dQdx_elife/xcorr;
}
