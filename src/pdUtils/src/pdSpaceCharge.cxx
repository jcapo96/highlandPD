#define pdSpaceCharge_C

#include "pdSpaceCharge.hxx"
#include "Parameters.hxx"

//********************************************************************
pdspacecharge::pdSpaceCharge::pdSpaceCharge(){
//********************************************************************

  //get parameters from parameters file
  fEnableSimSpatialSCE = ND::params().GetParameterI("pdUtils.SCE.EnableSimSpatialSCE");
  fEnableSimEfieldSCE  = ND::params().GetParameterI("pdUtils.SCE.EnableSimEfieldSCE");
  fEnableCalSpatialSCE = ND::params().GetParameterI("pdUtils.SCE.EnableCalSpatialSCE");
  fEnableCalEfieldSCE  = ND::params().GetParameterI("pdUtils.SCE.EnableCalEfieldSCE");
  fEfield = ND::params().GetParameterI("pdUtils.SCE.Efield");
  
  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true)){
    
    //read input file
    std::string fInputFilename = std::string(getenv("PIONANALYSISROOT"))+"/data/SCE_DataDriven_180kV_v3.root";
    infile = new TFile(fInputFilename.c_str(), "OPEN" );

    //check it exits
    if(!infile){
      std::cout << "Can't open file " << fInputFilename << std::endl;
      exit(1);
    }

    //Load in files
    hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
    hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
    hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
    hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
    hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
    hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
    
    hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
    hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
    hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
    hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
    hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
    hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
    
    hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
    hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
    hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
    hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
    hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
    hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");
    
    hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
    hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
    hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
    hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
    hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
    hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");
      
    hDx_sim_pos->SetDirectory(0);
    hDy_sim_pos->SetDirectory(0);
    hDz_sim_pos->SetDirectory(0);
    hEx_sim_pos->SetDirectory(0);
    hEy_sim_pos->SetDirectory(0);
    hEz_sim_pos->SetDirectory(0);
    
    hDx_sim_neg->SetDirectory(0);
    hDy_sim_neg->SetDirectory(0);
    hDz_sim_neg->SetDirectory(0);
    hEx_sim_neg->SetDirectory(0);
    hEy_sim_neg->SetDirectory(0);
    hEz_sim_neg->SetDirectory(0);
      
    SCEhistograms = {hDx_sim_pos, hDy_sim_pos, hDz_sim_pos, hEx_sim_pos, hEy_sim_pos, hEz_sim_pos, hDx_sim_neg, hDy_sim_neg, hDz_sim_neg, hEx_sim_neg, hEy_sim_neg, hEz_sim_neg};
  }
  
  if((fEnableCalSpatialSCE == true) || (fEnableCalEfieldSCE == true)){

    //read input file
    std::string fCalInputFilename = std::string(getenv("PIONANALYSISROOT"))+"/data/SCE_DataDriven_180kV_v3.root";
    infile = new TFile(fCalInputFilename.c_str(), "OPEN" );

    //check it exists
    if(!infile){
      std::cout << "Can't open file " << fCalInputFilename << std::endl;
      exit(1);
    }
    
    //Load in files
    hDx_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_X_Pos");
    hDy_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Y_Pos");
    hDz_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Z_Pos");
    hEx_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
    hEy_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
    hEz_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
    
    hDx_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_X_Neg");
    hDy_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Y_Neg");
    hDz_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Z_Neg");
    hEx_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
    hEy_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
    hEz_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
    
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
    
    CalSCEhistograms = {hDx_cal_pos, hDy_cal_pos, hDz_cal_pos, hEx_cal_pos, hEy_cal_pos, hEz_cal_pos, hDx_cal_neg, hDy_cal_neg, hDz_cal_neg, hEx_cal_neg, hEy_cal_neg, hEz_cal_neg};
  }
  infile->Close();
}

//********************************************************************
pdspacecharge::pdSpaceCharge::~pdSpaceCharge(){
//********************************************************************

  if(infile)infile->Delete();
}

//********************************************************************
bool pdspacecharge::pdSpaceCharge::EnableSimSpatialSCE() const{
//********************************************************************
   return fEnableSimSpatialSCE;
 }

//********************************************************************
bool pdspacecharge::pdSpaceCharge::EnableSimEfieldSCE() const{
//********************************************************************
  return fEnableSimEfieldSCE;
}

//********************************************************************
bool pdspacecharge::pdSpaceCharge::EnableCalSpatialSCE() const{
//********************************************************************
  return fEnableCalSpatialSCE;
}
 
//********************************************************************
bool pdspacecharge::pdSpaceCharge::EnableCalEfieldSCE() const{
//********************************************************************
  return fEnableCalEfieldSCE;
}

//********************************************************************
std::vector<double> pdspacecharge::pdSpaceCharge::GetPosOffsets(std::vector<double> const& tmp_point) const{
//********************************************************************  
  std::vector<double> thePosOffsets;
  std::vector<double> point = tmp_point;
  if(IsTooFarFromBoundaries(point)){
    thePosOffsets.resize(3,0.0);
    for(int i = 0; i < (int)thePosOffsets.size(); i++)point.push_back(-thePosOffsets[i]);
    return point;
  }
   if(!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point))point = PretendAtBoundary(point);
   
   if (point[0] > 0.){
     thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2));
     thePosOffsets[0] = thePosOffsets[0];
   } 
   else{
     thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(6), SCEhistograms.at(7), SCEhistograms.at(8));
     thePosOffsets[0] = -1.0*thePosOffsets[0];
   }
       
   return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//********************************************************************  
std::vector<double> pdspacecharge::pdSpaceCharge::GetCalPosOffsets(std::vector<double> const& tmp_point, int const& TPCid) const{
//********************************************************************  
  std::vector<double> thePosOffsets;
  std::vector<double> point = tmp_point;
   
  if(IsTooFarFromBoundaries(point)){
    thePosOffsets.resize(3,0.0);
    for(int i = 0; i < (int)thePosOffsets.size(); i++)point.push_back(-thePosOffsets[i]);
    return point;
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)){ 
    point = PretendAtBoundary(point); 
  }
  
  if ((TPCid == 2 || TPCid == 6 || TPCid == 10) && point[0]>-20.){
    if (point[0]<0.) point = {0.00001, point[1], point[2]};
    thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2));
    thePosOffsets[0] = -1.0*thePosOffsets[0];
  } 
  else if((TPCid == 1 || TPCid == 5 || TPCid == 9) && point[0]<20.) {
    if (point[0]>0.) point = {-0.00001, point[1], point[2]};
    thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(6), CalSCEhistograms.at(7), CalSCEhistograms.at(8));
    thePosOffsets[0] = -1.0*thePosOffsets[0];
  } 
  else thePosOffsets = {0., 0., 0.};
  
  return thePosOffsets;
}


//********************************************************************
std::vector<double> pdspacecharge::pdSpaceCharge::GetOffsetsVoxel(std::vector<double> const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const{
//********************************************************************
  return {hX->Interpolate(point[0],point[1],point[2]),
          hY->Interpolate(point[0],point[1],point[2]),
          hZ->Interpolate(point[0],point[1],point[2])};
}

//********************************************************************
std::vector<double> pdspacecharge::pdSpaceCharge::GetEfieldOffsets(std::vector<double> const& tmp_point) const{
//********************************************************************
  std::vector<double> theEfieldOffsets;
  std::vector<double> point = tmp_point;

  if(IsTooFarFromBoundaries(point)){
    theEfieldOffsets.resize(3,0.0);
    for(int i = 0; i < (int)theEfieldOffsets.size(); i++)point.push_back(-theEfieldOffsets[i]);
    return point;
  }
  if(!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if (point[0] > 0.) theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5));
  else               theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(9), SCEhistograms.at(10), SCEhistograms.at(11));
    
  return theEfieldOffsets;
}

//********************************************************************
std::vector<double> pdspacecharge::pdSpaceCharge::GetCalEfieldOffsets(std::vector<double> const& tmp_point, int const& TPCid) const{
//********************************************************************
  std::vector<double> theEfieldOffsets;
  std::vector<double> point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    theEfieldOffsets.resize(3,0.0);
    for(int i = 0; i < (int)theEfieldOffsets.size(); i++)point.push_back(-theEfieldOffsets[i]);
    return point;
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if ((TPCid == 2 || TPCid == 6 || TPCid == 10) && point[0] > -20.){
    if (point[0] < 0.) point = {0.00001, point[1], point[2]};
    theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5));
  }
  else if ((TPCid == 1 || TPCid == 5 || TPCid == 9) && point[0] < 20.){
    if (point[0] > 0.) point = {-0.00001, point[1], point[2]};
    theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(9), CalSCEhistograms.at(10), CalSCEhistograms.at(11));
  } 
  else theEfieldOffsets = {0., 0., 0.};
  
  return theEfieldOffsets;
}

//********************************************************************
double pdspacecharge::pdSpaceCharge::TransformX(double xVal) const{
//********************************************************************
  double xValNew;
  xValNew = (fabs(xVal)/100.0);
  //xValNew -= 1.8;
  return xValNew;
}

//********************************************************************
double pdspacecharge::pdSpaceCharge::TransformY(double yVal) const{
//********************************************************************
  double yValNew;
  yValNew = (6.00/6.08)*((yVal+0.2)/100.0);
  //yValNew -= 3.0;
  return yValNew;
}

//********************************************************************
double pdspacecharge::pdSpaceCharge::TransformZ(double zVal) const{
//********************************************************************
  double zValNew;
  zValNew = (7.20/6.97)*((zVal+0.8)/100.0);
  return zValNew;
}

//********************************************************************
bool pdspacecharge::pdSpaceCharge::IsInsideBoundaries(std::vector<double> const& point) const{
//********************************************************************
  return !((TMath::Abs(point[0]) <= 0.0)  || (TMath::Abs(point[0]) >= 360.0)
	   || (point[1]          <= 5.2)  || (point[1]             >= 604.0)
	   || (point[2]          <= -0.5) || (point[2]             >= 695.3)
	   );
} 

//********************************************************************
bool pdspacecharge::pdSpaceCharge::IsTooFarFromBoundaries(std::vector<double> const& point) const{
//********************************************************************  
  return((TMath::Abs(point[0]) < -20.0) || (TMath::Abs(point[0])  >= 360.0)
	 || (point[1]          < -14.8) || (point[1]              >  624.0)
	 || (point[2]          < -20.5) || (point[2]              >  715.3)
	 );
}
 
//********************************************************************
std::vector<double> pdspacecharge::pdSpaceCharge::PretendAtBoundary(std::vector<double> const& point) const{
//********************************************************************
  std::vector<double> rpoint;
  for(int i = 0; i < (int)point.size(); i++)rpoint.push_back(point[i]);
   
  if      (TMath::Abs(point[0]) ==    0.0    ) rpoint[0] = -0.00001;
  else if (TMath::Abs(point[0]) <     0.00001) rpoint[0] = TMath::Sign(point[0],1)*0.00001; 
  else if (TMath::Abs(point[0]) >=    360.0  ) rpoint[0] = TMath::Sign(point[0],1)*359.99999;
  
  if      (point[1] <=   5.2)                  rpoint[1] =   5.20001;
  else if (point[1] >= 604.0)                  rpoint[1] = 603.99999;
  
  if      (point[2] <=   -0.5)                 rpoint[2] =   -0.49999;
  else if (point[2] >= 695.3)                  rpoint[2] = 695.29999;
  
  return rpoint;
}
