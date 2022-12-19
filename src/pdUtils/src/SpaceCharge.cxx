#define SpaceCharge_C

#include "SpaceCharge.hxx"
#include "Parameters.hxx"

//********************************************************************
SpaceCharge::SpaceCharge(){
//********************************************************************

  fEnableElectronDiverterDistortions1 = false;
  fEnableElectronDiverterDistortions2 = false;
  fEDZCenterMin                       = -999.;
  fEDZCenterMax                       = -999.;
  fEDAXPosOffs                        = -999.;
  fEDBZPosOffs                        = -999.;
  fEDs                                = -999.;
  fEDChargeLossZLowMin                = -999.;
  fEDChargeLossZLowMax                = -999.;
  fEDChargeLossZHighMin               = -999.;
  fEDChargeLossZHighMax               = -999.;

  hDx_cal_pos = NULL;
  hDy_cal_pos = NULL;
  hDz_cal_pos = NULL;
  hEx_cal_pos = NULL;
  hEy_cal_pos = NULL;
  hEz_cal_pos = NULL;
  	   
  hDx_cal_neg = NULL;
  hDy_cal_neg = NULL;
  hDz_cal_neg = NULL;
  hEx_cal_neg = NULL;
  hEy_cal_neg = NULL;
  hEz_cal_neg = NULL;

  nominal_hDx_cal_pos = NULL;
  nominal_hDy_cal_pos = NULL;
  nominal_hDz_cal_pos = NULL;
  nominal_hEx_cal_pos = NULL;
  nominal_hEy_cal_pos = NULL;
  nominal_hEz_cal_pos = NULL;

  nominal_hDx_cal_neg = NULL;
  nominal_hDy_cal_neg = NULL;
  nominal_hDz_cal_neg = NULL;
  nominal_hEx_cal_neg = NULL;
  nominal_hEy_cal_neg = NULL;
  nominal_hEz_cal_neg = NULL;

  spline_dx_fwd_neg.clear();
  spline_dy_fwd_neg.clear();
  spline_dz_fwd_neg.clear();
  
  spline_dx_bkwd_neg.clear();
  spline_dy_bkwd_neg.clear();
  spline_dz_bkwd_neg.clear();
  
  spline_dEx_neg.clear();
  spline_dEy_neg.clear();
  spline_dEz_neg.clear();
  
  spline_dx_fwd_pos.clear();
  spline_dy_fwd_pos.clear();
  spline_dz_fwd_pos.clear();
  
  spline_dx_bkwd_pos.clear();
  spline_dy_bkwd_pos.clear();
  spline_dz_bkwd_pos.clear();
  
  spline_dEx_pos.clear();
  spline_dEy_pos.clear();
  spline_dEz_pos.clear();
  
  SCEhistograms.clear();

  _IsVaried = false;
}

//********************************************************************
SpaceCharge::~SpaceCharge(){
//********************************************************************

  delete hDx_cal_pos;
  delete hDy_cal_pos;
  delete hDz_cal_pos;
  delete hEx_cal_pos;
  delete hEy_cal_pos;
  delete hEz_cal_pos;

  delete hDx_cal_neg;
  delete hDy_cal_neg;
  delete hDz_cal_neg;
  delete hEx_cal_neg;
  delete hEy_cal_neg;
  delete hEz_cal_neg;

  delete nominal_hDx_cal_pos;
  delete nominal_hDy_cal_pos;
  delete nominal_hDz_cal_pos;
  delete nominal_hEx_cal_pos;
  delete nominal_hEy_cal_pos;
  delete nominal_hEz_cal_pos;

  delete nominal_hDx_cal_neg;
  delete nominal_hDy_cal_neg;
  delete nominal_hDz_cal_neg;
  delete nominal_hEx_cal_neg;
  delete nominal_hEy_cal_neg;
  delete nominal_hEz_cal_neg;

  ClearSplines();

  for(int i = 0; i < (int)SCEhistograms.size(); i++)
    delete SCEhistograms[i];
  SCEhistograms.clear();
}

//********************************************************************
void SpaceCharge::Initialize(){
//********************************************************************

  fEnableElectronDiverterDistortions1 = true;//ND::params().GetParameterF("EnableElectronDiverterDistortions1");
  fEnableElectronDiverterDistortions2 = true;//ND::params().GetParameterF("EnableElectronDiverterDistortions2");
  fEDZCenterMin                       = 231.;//ND::params().GetParameterI("EDZCenter");
  fEDZCenterMax                       = 463.;//ND::params().GetParameterI("EDZCenter");
  fEDAXPosOffs                        = -3.;//ND::params().GetParameterF("EDAXPosOffs");
  fEDBZPosOffs                        = -1;//ND::params().GetParameterF("EDBZPosOffs");
  fEDs                                = 2.5;//ND::params().GetParameterF("EDs");
  fEDChargeLossZLowMin                = 229.;//ND::params().GetParameterF("EDChargeLossZLowMin");
  fEDChargeLossZLowMax                = 461.;//ND::params().GetParameterF("EDChargeLossZLowMax");
  fEDChargeLossZHighMin               = 233.;//ND::params().GetParameterF("EDChargeLossZHighMin");
  fEDChargeLossZHighMax               = 465.;//ND::params().GetParameterF("EDChargeLossZHighMax");
  
  fInputFilename = std::string(getenv("PDUTILSROOT"))+"/data/SCE_DataDriven_180kV_v4.root";
  fInputFile = new TFile(fInputFilename.c_str(), "OPEN");
  
  //check it exits
  if(!fInputFile){
    std::cout << "Can't open file " << fInputFilename << std::endl;
    exit(1);
  }

  InitializeHistograms();
  
  InitializeSplines();
    
  CalSCEhistograms = {hDx_cal_pos, hDy_cal_pos, hDz_cal_pos, hEx_cal_pos, hEy_cal_pos, hEz_cal_pos, hDx_cal_neg, hDy_cal_neg, hDz_cal_neg, hEx_cal_neg, hEy_cal_neg, hEz_cal_neg};

  fInputFile->Close();
}


//********************************************************************
TVector3 SpaceCharge::GetPosOffsets(TVector3 const& point) const {
//********************************************************************
  
  TVector3 Offsets(0,0,0);
  TVector3 p = point;
 
  if(IsTooFarFromBoundaries(p))return Offsets;

  if(!IsInsideBoundaries(p)&&!IsTooFarFromBoundaries(p)) p = PretendAtBoundary(p);
  
  if(p.X() > 0.){
    Offsets = GetOffsets(p, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2), 1, 2);
    Offsets.SetX(-1.0*Offsets.X());
  } 
  else{
    Offsets = GetOffsets(p, SCEhistograms.at(6), SCEhistograms.at(7), SCEhistograms.at(8), 1, 1);
    Offsets.SetX(-1.0*Offsets.X());
  }
       
  TVector3 pafteroffset(p.X()+Offsets.X(), p.Y()+Offsets.Y(), p.Z()+Offsets.Z());
  TVector3 edoffset = ElectronDiverterPosOffsets(pafteroffset);
  Offsets.SetX(Offsets.X()+edoffset.X());
  Offsets.SetY(Offsets.Y()+edoffset.Y());
  Offsets.SetZ(Offsets.Z()+edoffset.Z());

  return Offsets;
}

//********************************************************************
TVector3 SpaceCharge::GetCalPosOffsets(TVector3 const& point, int const& TPCid) const {
//********************************************************************

  TVector3 Offsets(0,0,0);
  TVector3 p = point;

  if(IsTooFarFromBoundaries(p)) return Offsets;

  if(!IsInsideBoundaries(p) && !IsTooFarFromBoundaries(p))
  	p = PretendAtBoundary(p); 
  
  if((TPCid == 2 || TPCid == 6 || TPCid == 10) && p.X()>-20.){
    if(p.X()<0.) p = {0.00001, p.Y(), p.Z()};
    Offsets = GetOffsets(p, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2), 2, 2);
    Offsets[0] = -1.0*Offsets[0];
  } 
  else if((TPCid == 1 || TPCid == 5 || TPCid == 9) && p.X()<20.) {
    if(p.X()>0.) p= {-0.00001, p.Y(), p.Z()};
    Offsets = GetOffsets(p, CalSCEhistograms.at(6), CalSCEhistograms.at(7), CalSCEhistograms.at(8), 2, 1);
    Offsets[0] = -1.0*Offsets[0];
  } 

  return Offsets;
}

//********************************************************************
TVector3 SpaceCharge::GetOffsets(TVector3 const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const {
//********************************************************************

  return {InterpolateSplines(hX,point.X(),point.Y(),point.Z(),1,maptype,driftvol),
          InterpolateSplines(hY,point.X(),point.Y(),point.Z(),2,maptype,driftvol),
          InterpolateSplines(hZ,point.X(),point.Y(),point.Z(),3,maptype,driftvol)};
}

//********************************************************************
TVector3 SpaceCharge::GetEfieldOffsets(TVector3 const& point) const {
//********************************************************************
  
  TVector3 Offsets(0,0,0);
  TVector3 p = point;

  if(IsTooFarFromBoundaries(p))return Offsets;

  if(!IsInsideBoundaries(p) && !IsTooFarFromBoundaries(p)) p = PretendAtBoundary(p);
  
  if(p.X() > 0.)Offsets = GetOffsets(p, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5), 3, 2);
  else          Offsets = GetOffsets(p, SCEhistograms.at(9), SCEhistograms.at(10),SCEhistograms.at(11),3, 1);
  Offsets.SetX(-1.0*Offsets.X());
  Offsets.SetY(-1.0*Offsets.Y());
  Offsets.SetZ(-1.0*Offsets.Z());
    
  return Offsets;
}

//********************************************************************
TVector3 SpaceCharge::GetCalEfieldOffsets(TVector3 const& point, int const& TPCid) const { 
//********************************************************************
 
 TVector3 Offsets(0,0,0);
 TVector3 p = point;

 if(IsTooFarFromBoundaries(p))return Offsets;
 
 if(!IsInsideBoundaries(p) && !IsTooFarFromBoundaries(p)) p = PretendAtBoundary(p);
  
 if((TPCid == 2 || TPCid == 6 || TPCid == 10) && p.X()>-20.){
   if(p.X()<0.) p = {0.00001, p.Y(), p.Z()};
   Offsets= GetOffsets(p, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5), 3, 2);
 }
 else if((TPCid == 1 || TPCid == 5 || TPCid == 9) && p.X()<20.){
   if(p.X()>0.) p = {-0.00001, p.Y(), p.Z()};
   Offsets = GetOffsets(p, CalSCEhistograms.at(9), CalSCEhistograms.at(10), CalSCEhistograms.at(11), 3, 1);
 } 
 Offsets.SetX(-1.0*Offsets.X());
 Offsets.SetY(-1.0*Offsets.Y());
 Offsets.SetZ(-1.0*Offsets.Z());
 
 return Offsets;
}


//********************************************************************
bool SpaceCharge::IsInsideBoundaries(TVector3 const& point) const {
//********************************************************************
  return !(
	   (TMath::Abs(point.X()) <= 0.0) || (TMath::Abs(point.X()) >= 360.0)
	   || (point.Y()             <= 5.2) || (point.Y()             >= 604.0)
	   || (point.Z()             <= -0.5) || (point.Z()             >= 695.3)
	   );
} 
  
//********************************************************************
bool SpaceCharge::IsTooFarFromBoundaries(TVector3 const& point) const {
//********************************************************************
  
  return (
	  (TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X())  >= 360.0)
	  || (point.Y()             < -14.8) || (point.Y()              >  624.0)
	  || (point.Z()             < -20.5) || (point.Z()              >  715.3)
	  );
}

//********************************************************************
TVector3 SpaceCharge::PretendAtBoundary(TVector3 const& point) const {
//********************************************************************
  
  double x = point.X(), y = point.Y(), z = point.Z();
  
  if      (TMath::Abs(point.X()) ==    0.0    )   x = -0.00001;
  else if (TMath::Abs(point.X()) <	 0.00001) x = TMath::Sign(point.X(),1)*0.00001; 
  else if (TMath::Abs(point.X()) >=    360.0  )   x = TMath::Sign(point.X(),1)*359.99999;
  
  if      (point.Y() <=   5.2) y = 5.20001;
  else if (point.Y() >= 604.0) y = 603.99999;
  
  if      (point.Z() <=   -0.5) z = -0.49999;
  else if (point.Z() >= 695.3)  z = 695.29999;
   
  return {x, y, z};
}

//********************************************************************
TSpline3* SpaceCharge::MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const {
//********************************************************************

  TSpline3 *spline = 0;
  
  std::vector<double> a, b;
  if (dim1 == 1) {
    for (int x = 1; x <= spline_hist->GetNbinsX(); ++x) {
      a.push_back(spline_hist->GetXaxis()->GetBinCenter(x));
      b.push_back(spline_hist->GetBinContent(x, dim2_bin, dim3_bin));
    }
  }
  else if (dim1 == 2) {
    for(int y = 1; y <= spline_hist->GetNbinsY(); ++y) {
      a.push_back(spline_hist->GetYaxis()->GetBinCenter(y));
      b.push_back(spline_hist->GetBinContent(dim2_bin, y, dim3_bin));
    }
  }
  else if (dim1 == 3) {
    for (int z = 1; z <= spline_hist->GetNbinsZ(); z++) {
      a.push_back(spline_hist->GetZaxis()->GetBinCenter(z));
      b.push_back(spline_hist->GetBinContent(dim2_bin, dim3_bin, z));
    }
  }
  else {
    std::cout << "Unkown dimension " << dim1 << std::endl;
    std::exit(1);
  }

  if(maptype == 1)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }
  else if(maptype == 2)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }
  else if(maptype == 3)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }


  return spline;
}

//********************************************************************
double SpaceCharge::InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const {
//********************************************************************

  //std::cout << "Interpolating " << interp_hist->GetName() << std::endl;
  int bin_x = interp_hist->GetXaxis()->FindBin(xVal);
  int bin_y = interp_hist->GetYaxis()->FindBin(yVal);
  int bin_z = interp_hist->GetZaxis()->FindBin(zVal);

  int bincenter_x = interp_hist->GetXaxis()->GetBinCenter(bin_x);
  int bincenter_y = interp_hist->GetYaxis()->GetBinCenter(bin_y);
  int bincenter_z = interp_hist->GetZaxis()->GetBinCenter(bin_z);

  int max_x = interp_hist->GetNbinsX();
  int max_y = interp_hist->GetNbinsY();
  int max_z = interp_hist->GetNbinsZ();
  
  int low_x;
  int high_x;
  if(bin_x <= 1)
  {
    low_x = 1;
    high_x = 2;
  }
  else if(bin_x >= max_x)
  {
    low_x = max_x-1;
    high_x = max_x;
  }
  else if(xVal > bincenter_x)
  {
    low_x = bin_x;
    high_x = bin_x+1;
  }
  else
  {
    low_x = bin_x-1;
    high_x = bin_x;
  }

  int low_y;
  int high_y;
  if(bin_y <= 1)
  {
    low_y = 1;
    high_y = 2;
  }
  else if(bin_y >= max_y)
  {
    low_y = max_y-1;
    high_y = max_y;
  }
  else if(yVal > bincenter_y)
  {
    low_y = bin_y;
    high_y = bin_y+1;
  }
  else
  {
    low_y = bin_y-1;
    high_y = bin_y;
  }

  int low_z;
  int high_z;
  if(bin_z <= 1)
  {
    low_z = 1;
    high_z = 2;
  }
  else if(bin_z >= max_z)
  {
    low_z = max_z-1;
    high_z = max_z;
  }
  else if(zVal > bincenter_z)
  {
    low_z = bin_z;
    high_z = bin_z+1;
  }
  else
  {
    low_z = bin_z-1;
    high_z = bin_z;
  }

  double interp_val = 0.0;
  
  if(dim == 1)
  {
    double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_neg[high_y-1][high_z-1]->Eval(xVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_pos[high_y-1][high_z-1]->Eval(xVal);
      }
    }

    interp_val = (f_11*(a_2-yVal)*(b_2-zVal) + f_21*(yVal-a_1)*(b_2-zVal) + f_12*(a_2-yVal)*(zVal-b_1) + f_22*(yVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 2)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_neg[high_x-1][high_z-1]->Eval(yVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_pos[high_x-1][high_z-1]->Eval(yVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-zVal) + f_21*(xVal-a_1)*(b_2-zVal) + f_12*(a_2-xVal)*(zVal-b_1) + f_22*(xVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 3)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_neg[high_x-1][high_y-1]->Eval(zVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_pos[high_x-1][high_y-1]->Eval(zVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-yVal) + f_21*(xVal-a_1)*(b_2-yVal) + f_12*(a_2-xVal)*(yVal-b_1) + f_22*(xVal-a_1)*(yVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }

  return interp_val;
}

//********************************************************************
TVector3 SpaceCharge::ElectronDiverterPosOffsets(TVector3 const& point) const {
//********************************************************************

  double z = point.Z();
  TVector3 Offsets(0,0,0);

  if (!fEnableElectronDiverterDistortions1){
    if(point.X()<0){
      if(z > fEDChargeLossZLowMin && z < fEDChargeLossZHighMin)
	Offsets.SetXYZ(2E9,2E9,2E9);
      else{
	double zdiff = z - fEDZCenterMin;
	double zexp = TMath::Exp( -TMath::Sq(zdiff/fEDs) );
	Offsets.SetZ(Offsets.Z() + fEDBZPosOffs * zdiff * zexp);
	
	// the timing offsets need to be computed after the z shift
	double zdiffc = zdiff + Offsets.Z();
	double zexpc = TMath::Exp( -TMath::Sq(zdiffc/fEDs) );
	Offsets.SetX(Offsets.X() + fEDAXPosOffs * zexpc);
      }
    }
  }

  if (!fEnableElectronDiverterDistortions2){
    if(point.X()<0){
      if(z > fEDChargeLossZLowMax && z < fEDChargeLossZHighMax)
	Offsets.SetXYZ(2E9,2E9,2E9);
      else{
	double zdiff = z - fEDZCenterMax;
	double zexp = TMath::Exp( -TMath::Sq(zdiff/fEDs) );
	Offsets.SetZ(Offsets.Z() + fEDBZPosOffs * zdiff * zexp);
	
	// the timing offsets need to be computed after the z shift
	double zdiffc = zdiff + Offsets.Z();
	double zexpc = TMath::Exp( -TMath::Sq(zdiffc/fEDs) );
	Offsets.SetX(Offsets.X() + fEDAXPosOffs * zexpc);
      }
    }
  }
  return Offsets;
}

//********************************************************************
void SpaceCharge::ApplyPositionCorrection(AnaParticlePD* part) const {
//********************************************************************

  if(!part)return;
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    ApplyPositionCorrection(part->Hits[2][ihit]);
}

//********************************************************************
void SpaceCharge::ApplyPositionCorrection(AnaHitPD& hit) const {
//********************************************************************

  TVector3 offset = GetCalPosOffsets(hit.Position_NoSCE, hit.TPCid); //*-1 check this
  hit.Position.SetX(hit.Position_NoSCE.X() - offset.X());
  hit.Position.SetY(hit.Position_NoSCE.Y() + offset.Y());
  hit.Position.SetZ(hit.Position_NoSCE.Z() + offset.Z());
}

//********************************************************************
void SpaceCharge::ApplyDisplacementVariation(const double var){
//********************************************************************
    
  int nBinsX = hDx_cal_pos->GetNbinsX();
  int nBinsY = hDx_cal_pos->GetNbinsY();
  int nBinsZ = hDx_cal_pos->GetNbinsZ();
  
  //loop over bins of the histograms
  for(int x = 1; x <= nBinsX; x++){
    for(int y = 1; y <= nBinsY; y++){
      for(int z = 1; z <= nBinsZ; z++){
  	hDx_cal_pos->SetBinContent(x,y,z,nominal_hDx_cal_pos->GetBinContent(x,y,z)*var);
  	hDy_cal_pos->SetBinContent(x,y,z,nominal_hDy_cal_pos->GetBinContent(x,y,z)*var);
  	hDz_cal_pos->SetBinContent(x,y,z,nominal_hDz_cal_pos->GetBinContent(x,y,z)*var);

  	hDx_cal_neg->SetBinContent(x,y,z,nominal_hDx_cal_neg->GetBinContent(x,y,z)*var);
  	hDy_cal_neg->SetBinContent(x,y,z,nominal_hDy_cal_neg->GetBinContent(x,y,z)*var);
  	hDz_cal_neg->SetBinContent(x,y,z,nominal_hDz_cal_neg->GetBinContent(x,y,z)*var);
      }
    }
  }

  // hDx_cal_pos->Scale(var);
  // hDy_cal_pos->Scale(var);
  // hDz_cal_pos->Scale(var);

  // hDx_cal_neg->Scale(var);
  // hDy_cal_neg->Scale(var);
  // hDz_cal_neg->Scale(var);
    
  //prepare new splines
  ResetSplines();

  //set true
  _IsVaried = true;
}

//********************************************************************
void SpaceCharge::ResetToNominal(){
//********************************************************************
    
  delete hDx_cal_pos;
  delete hDy_cal_pos;
  delete hDz_cal_pos;
  delete hEx_cal_pos;
  delete hEy_cal_pos;
  delete hEz_cal_pos;

  hDx_cal_pos = (TH3F*)nominal_hDx_cal_pos->Clone();
  hDy_cal_pos = (TH3F*)nominal_hDy_cal_pos->Clone();
  hDz_cal_pos = (TH3F*)nominal_hDz_cal_pos->Clone();
  hEx_cal_pos = (TH3F*)nominal_hEx_cal_pos->Clone();
  hEy_cal_pos = (TH3F*)nominal_hEy_cal_pos->Clone();
  hEz_cal_pos = (TH3F*)nominal_hEz_cal_pos->Clone();

  // int nBinsX = hDx_cal_pos->GetNbinsX();
  // int nBinsY = hDx_cal_pos->GetNbinsY();
  // int nBinsZ = hDx_cal_pos->GetNbinsZ();

  // //loop over bins of the histograms //this is not efficient but it does not crash
  // for(int x = 1; x <= nBinsX; x++){
  //   for(int y = 1; y <= nBinsY; y++){
  //     for(int z = 1; z <= nBinsZ; z++){
  // 	hDx_cal_pos->SetBinContent(x,y,z,nominal_hDx_cal_pos->GetBinContent(x,y,z));
  // 	hDy_cal_pos->SetBinContent(x,y,z,nominal_hDy_cal_pos->GetBinContent(x,y,z));
  // 	hDz_cal_pos->SetBinContent(x,y,z,nominal_hDz_cal_pos->GetBinContent(x,y,z));

  // 	hDx_cal_neg->SetBinContent(x,y,z,nominal_hDx_cal_neg->GetBinContent(x,y,z));
  // 	hDy_cal_neg->SetBinContent(x,y,z,nominal_hDy_cal_neg->GetBinContent(x,y,z));
  // 	hDz_cal_neg->SetBinContent(x,y,z,nominal_hDz_cal_neg->GetBinContent(x,y,z));
  //     }
  //   }
  // }

  //reset splines
  ResetSplines();
}

//********************************************************************
void SpaceCharge::ResetSplines(){
//********************************************************************
    
  ClearSplines();
  InitializeSplines();
}

//********************************************************************
void SpaceCharge::ClearSplines(){
//********************************************************************

  ClearVectorOfSplines(spline_dx_fwd_neg);
  ClearVectorOfSplines(spline_dy_fwd_neg);
  ClearVectorOfSplines(spline_dz_fwd_neg);

  ClearVectorOfSplines(spline_dx_bkwd_neg);
  ClearVectorOfSplines(spline_dy_bkwd_neg);
  ClearVectorOfSplines(spline_dz_bkwd_neg);

  ClearVectorOfSplines(spline_dEx_neg);
  ClearVectorOfSplines(spline_dEy_neg);
  ClearVectorOfSplines(spline_dEz_neg);

  ClearVectorOfSplines(spline_dx_fwd_pos);
  ClearVectorOfSplines(spline_dy_fwd_pos);
  ClearVectorOfSplines(spline_dz_fwd_pos);

  ClearVectorOfSplines(spline_dx_bkwd_pos);
  ClearVectorOfSplines(spline_dy_bkwd_pos);
  ClearVectorOfSplines(spline_dz_bkwd_pos);

  ClearVectorOfSplines(spline_dEx_pos);
  ClearVectorOfSplines(spline_dEy_pos);
  ClearVectorOfSplines(spline_dEz_pos);
}

//********************************************************************
void SpaceCharge::ClearVectorOfSplines(std::vector<std::vector<TSpline3*>> &spline){
//********************************************************************
  
  for(int i = 0; i < (int)spline.size(); i++){
    for(int j = 0; j < (int)spline[i].size(); j++)
      delete spline[i][j];
    spline[i].clear();
  }
  spline.clear();
}

//********************************************************************
void SpaceCharge::InitializeHistograms(){
//********************************************************************
  
  //Load in files
  TH3F* hDx_cal_pos_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_X_Pos");
  TH3F* hDy_cal_pos_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_Y_Pos");
  TH3F* hDz_cal_pos_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_Z_Pos");
  TH3F* hEx_cal_pos_orig = (TH3F*)fInputFile->Get("Reco_ElecField_X_Pos");
  TH3F* hEy_cal_pos_orig = (TH3F*)fInputFile->Get("Reco_ElecField_Y_Pos");
  TH3F* hEz_cal_pos_orig = (TH3F*)fInputFile->Get("Reco_ElecField_Z_Pos");
        
  TH3F* hDx_cal_neg_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_X_Neg");
  TH3F* hDy_cal_neg_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_Y_Neg");
  TH3F* hDz_cal_neg_orig = (TH3F*)fInputFile->Get("RecoBkwd_Displacement_Z_Neg");
  TH3F* hEx_cal_neg_orig = (TH3F*)fInputFile->Get("Reco_ElecField_X_Neg");
  TH3F* hEy_cal_neg_orig = (TH3F*)fInputFile->Get("Reco_ElecField_Y_Neg");
  TH3F* hEz_cal_neg_orig = (TH3F*)fInputFile->Get("Reco_ElecField_Z_Neg");
  
  hDx_cal_pos = (TH3F*)hDx_cal_pos_orig->Clone("hDx_pos");
  hDy_cal_pos = (TH3F*)hDy_cal_pos_orig->Clone("hDy_pos");
  hDz_cal_pos = (TH3F*)hDz_cal_pos_orig->Clone("hDz_pos");
  hEx_cal_pos = (TH3F*)hEx_cal_pos_orig->Clone("hEx_pos");
  hEy_cal_pos = (TH3F*)hEy_cal_pos_orig->Clone("hEy_pos");
  hEz_cal_pos = (TH3F*)hEz_cal_pos_orig->Clone("hEz_pos");
  
  hDx_cal_neg = (TH3F*)hDx_cal_neg_orig->Clone("hDx_neg");
  hDy_cal_neg = (TH3F*)hDy_cal_neg_orig->Clone("hDy_neg");
  hDz_cal_neg = (TH3F*)hDz_cal_neg_orig->Clone("hDz_neg");
  hEx_cal_neg = (TH3F*)hEx_cal_neg_orig->Clone("hEx_neg");
  hEy_cal_neg = (TH3F*)hEy_cal_neg_orig->Clone("hEy_neg");
  hEz_cal_neg = (TH3F*)hEz_cal_neg_orig->Clone("hEz_neg");

  nominal_hDx_cal_pos = (TH3F*)hDx_cal_pos_orig->Clone("hDx_pos");
  nominal_hDy_cal_pos = (TH3F*)hDy_cal_pos_orig->Clone("hDy_pos");
  nominal_hDz_cal_pos = (TH3F*)hDz_cal_pos_orig->Clone("hDz_pos");
  nominal_hEx_cal_pos = (TH3F*)hEx_cal_pos_orig->Clone("hEx_pos");
  nominal_hEy_cal_pos = (TH3F*)hEy_cal_pos_orig->Clone("hEy_pos");
  nominal_hEz_cal_pos = (TH3F*)hEz_cal_pos_orig->Clone("hEz_pos");

  nominal_hDx_cal_neg = (TH3F*)hDx_cal_neg_orig->Clone("hDx_neg");
  nominal_hDy_cal_neg = (TH3F*)hDy_cal_neg_orig->Clone("hDy_neg");
  nominal_hDz_cal_neg = (TH3F*)hDz_cal_neg_orig->Clone("hDz_neg");
  nominal_hEx_cal_neg = (TH3F*)hEx_cal_neg_orig->Clone("hEx_neg");
  nominal_hEy_cal_neg = (TH3F*)hEy_cal_neg_orig->Clone("hEy_neg");
  nominal_hEz_cal_neg = (TH3F*)hEz_cal_neg_orig->Clone("hEz_neg");

  hDx_cal_pos->SetDirectory(0);
  hDy_cal_pos->SetDirectory(0);
  hDz_cal_pos->SetDirectory(0);
  hEx_cal_pos->SetDirectory(0);
  hEy_cal_pos->SetDirectory(0);
  hEz_cal_pos->SetDirectory(0);
  
  hDx_cal_neg->SetDirectory(0);
  hDy_cal_neg->SetDirectory(0);
  hDz_cal_neg->SetDirectory(0);
  hEx_cal_neg->SetDirectory(0);
  hEy_cal_neg->SetDirectory(0);
  hEz_cal_neg->SetDirectory(0);
  
  nominal_hDx_cal_pos->SetDirectory(0);
  nominal_hDy_cal_pos->SetDirectory(0);
  nominal_hDz_cal_pos->SetDirectory(0);
  nominal_hEx_cal_pos->SetDirectory(0);
  nominal_hEy_cal_pos->SetDirectory(0);
  nominal_hEz_cal_pos->SetDirectory(0);

  nominal_hDx_cal_neg->SetDirectory(0);
  nominal_hDy_cal_neg->SetDirectory(0);
  nominal_hDz_cal_neg->SetDirectory(0);
  nominal_hEx_cal_neg->SetDirectory(0);
  nominal_hEy_cal_neg->SetDirectory(0);
  nominal_hEz_cal_neg->SetDirectory(0);
}

//********************************************************************
void SpaceCharge::InitializeSplines(){
//********************************************************************

  int nBinsX = hDx_cal_pos->GetNbinsX();
  int nBinsY = hDx_cal_pos->GetNbinsY();
  int nBinsZ = hDx_cal_pos->GetNbinsZ();
  
  for(int y = 1; y <= nBinsY; y++){
    spline_dx_bkwd_neg.push_back(std::vector<TSpline3*>());
    spline_dx_bkwd_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dx_bkwd_neg.back().push_back(MakeSpline(hDx_cal_neg,1,y,z,2,1));
      spline_dx_bkwd_pos.back().push_back(MakeSpline(hDx_cal_pos,1,y,z,2,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dy_bkwd_neg.push_back(std::vector<TSpline3*>());
    spline_dy_bkwd_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dy_bkwd_neg.back().push_back(MakeSpline(hDy_cal_neg,2,x,z,2,1));
      spline_dy_bkwd_pos.back().push_back(MakeSpline(hDy_cal_pos,2,x,z,2,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dz_bkwd_neg.push_back(std::vector<TSpline3*>());
    spline_dz_bkwd_pos.push_back(std::vector<TSpline3*>());
    for(int y = 1; y <= nBinsY; y++){
      spline_dz_bkwd_neg.back().push_back(MakeSpline(hDz_cal_neg,3,x,y,2,1));
      spline_dz_bkwd_pos.back().push_back(MakeSpline(hDz_cal_pos,3,x,y,2,2));
    }
  }
  nBinsX = hEx_cal_neg->GetNbinsX();
  nBinsY = hEx_cal_neg->GetNbinsY();
  nBinsZ = hEx_cal_neg->GetNbinsZ();
  for(int y = 1; y <= nBinsY; y++){
    spline_dEx_neg.push_back(std::vector<TSpline3*>());
    spline_dEx_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dEx_neg.back().push_back(MakeSpline(hEx_cal_neg,1,y,z,3,1));
      spline_dEx_pos.back().push_back(MakeSpline(hEx_cal_pos,1,y,z,3,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dEy_neg.push_back(std::vector<TSpline3*>());
    spline_dEy_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dEy_neg.back().push_back(MakeSpline(hEy_cal_neg,2,x,z,3,1));
      spline_dEy_pos.back().push_back(MakeSpline(hEy_cal_pos,2,x,z,3,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dEz_neg.push_back(std::vector<TSpline3*>());
    spline_dEz_pos.push_back(std::vector<TSpline3*>());
    for(int y = 1; y <= nBinsY; y++){
      spline_dEz_neg.back().push_back(MakeSpline(hEz_cal_neg,3,x,y,3,1));
      spline_dEz_pos.back().push_back(MakeSpline(hEz_cal_pos,3,x,y,3,2));
    }
  }

  /*//Load in files
  TH3F* hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
  TH3F* hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
  TH3F* hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
  TH3F* hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
  TH3F* hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
  TH3F* hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
    
  TH3F* hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
  TH3F* hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
  TH3F* hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
  TH3F* hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
  TH3F* hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
  TH3F* hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
  
  TH3F* hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
  TH3F* hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
  TH3F* hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
  TH3F* hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
  TH3F* hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
  TH3F* hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");
  
  TH3F* hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
  TH3F* hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
  TH3F* hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
  TH3F* hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
  TH3F* hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
  TH3F* hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");
  
  int nBinsX = hDx_sim_pos_orig->GetNbinsX();
  int nBinsY = hDx_sim_pos_orig->GetNbinsY();
  int nBinsZ = hDx_sim_pos_orig->GetNbinsZ();
  for(int y = 1; y <= nBinsY; y++){
    spline_dx_fwd_neg.push_back(std::vector<TSpline3*>());
    spline_dx_fwd_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dx_fwd_neg.back().push_back(MakeSpline(hDx_sim_neg,1,y,z,1,1));
      spline_dx_fwd_pos.back().push_back(MakeSpline(hDx_sim_pos,1,y,z,1,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dy_fwd_neg.push_back(std::vector<TSpline3*>());
    spline_dy_fwd_pos.push_back(std::vector<TSpline3*>());
    
    for(int z = 1; z <= nBinsZ; z++){
      spline_dy_fwd_neg.back().push_back(MakeSpline(hDy_sim_neg,2,x,z,1,1));
      spline_dy_fwd_pos.back().push_back(MakeSpline(hDy_sim_pos,2,x,z,1,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dz_fwd_neg.push_back(std::vector<TSpline3*>());
    spline_dz_fwd_pos.push_back(std::vector<TSpline3*>());
    for(int y = 1; y <= nBinsY; y++){
      spline_dz_fwd_neg.back().push_back(MakeSpline(hDz_sim_neg,3,x,y,1,1));
      spline_dz_fwd_pos.back().push_back(MakeSpline(hDz_sim_pos,3,x,y,1,2));
    }
  }
    
  nBinsX = hEx_sim_pos_orig->GetNbinsX();
  nBinsY = hEx_sim_pos_orig->GetNbinsY();
  nBinsZ = hEx_sim_pos_orig->GetNbinsZ();
  for(int y = 1; y <= nBinsY; y++){
    spline_dEx_neg.push_back(std::vector<TSpline3*>());
    spline_dEx_pos.push_back(std::vector<TSpline3*>());
    for(int z = 1; z <= nBinsZ; z++){
      spline_dEx_neg.back().push_back(MakeSpline(hEx_sim_neg,1,y,z,3,1));
      spline_dEx_pos.back().push_back(MakeSpline(hEx_sim_pos,1,y,z,3,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dEy_neg.push_back(std::vector<TSpline3*>());
    spline_dEy_pos.push_back(std::vector<TSpline3*>());
    
    for(int z = 1; z <= nBinsZ; z++){
      spline_dEy_neg.back().push_back(MakeSpline(hEy_sim_neg,2,x,z,3,1));
      spline_dEy_pos.back().push_back(MakeSpline(hEy_sim_pos,2,x,z,3,2));
    }
  }
  for(int x = 1; x <= nBinsX; x++){
    spline_dEz_neg.push_back(std::vector<TSpline3*>());
    spline_dEz_pos.push_back(std::vector<TSpline3*>());
    for(int y = 1; y <= nBinsY; y++){
      spline_dEz_neg.back().push_back(MakeSpline(hEz_sim_neg,3,x,y,3,1));
      spline_dEz_pos.back().push_back(MakeSpline(hEz_sim_pos,3,x,y,3,2));
    }
  }
  
  SCEhistograms = {hDx_sim_pos, hDy_sim_pos, hDz_sim_pos, hEx_sim_pos, hEy_sim_pos, hEz_sim_pos, hDx_sim_neg, hDy_sim_neg, hDz_sim_neg, hEx_sim_neg, hEy_sim_neg, hEz_sim_neg};
   
  infile->Close();  */
}

