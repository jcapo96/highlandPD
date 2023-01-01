//*********************************************************************
//minimal implementetation of the space charge effects in ProtoDUNE-SP
//only considering spline_th3 methods 
//and backwards (correction) implementation
//Get_Offsets give the effect forward of the SCE effect (it does not work) 
//GetCaL_Offsets give the effect backward (it works)
//
//M. García, migarpez@ific.uv.es
//*********************************************************************

#ifndef SpaceCharge_h
#define SpaceCharge_h

#include "TFile.h"
#include "TF1.h"
#include "TH3F.h"
#include "TSpline.h"
#include "TVector3.h"

#include "pdDataClasses.hxx"

class SpaceCharge{
public:
  
  SpaceCharge();
  virtual ~SpaceCharge();
    
  void Initialize();

  TVector3 GetPosOffsets(TVector3 const& point) const;
  TVector3 GetEfieldOffsets(TVector3 const& point) const;
  TVector3 GetCalPosOffsets(TVector3 const& point, int const& TPCid) const;
  TVector3 GetCalEfieldOffsets(TVector3 const& point, int const& TPCid) const;

  void ApplyPositionCorrection(AnaParticlePD* part) const;
  void ApplyPositionCorrection(AnaHitPD& hit) const;

  void ApplyDisplacementVariation(const double var);
  void ApplyVoxelVariation(UInt_t xbin, UInt_t ybin, UInt_t zbin, double var, bool reset_splines = true);
  void ResetToNominal();
  void ResetSplines();

  bool IsVaried() const {return _IsVaried;}

  UInt_t GetNbinsX() const;
  UInt_t GetNbinsY() const;
  UInt_t GetNbinsZ() const;

  double GetBinCenterX(UInt_t bin) const;
  double GetBinCenterY(UInt_t bin) const;
  double GetBinCenterZ(UInt_t bin) const;
  
protected:
  
  TVector3 GetOffsets(TVector3 const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const;

  void InitializeHistograms();
  void InitializeHistogram(TH3F** h, const char* histo, const char* newname);
  void InitializeSplines();

  void ClearSplines();
  void ClearVectorOfSplines(std::vector<std::vector<TSpline3*>> &splines);

  std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(12); //Histograms are Dx, Dy, Dz, dEx/E0, dEy/E0, dEz/E0 (positive; repeat for negative)
  std::vector<TH3F*> CalSCEhistograms = std::vector<TH3F*>(12); 
  
  bool IsInsideBoundaries(TVector3 const& point) const;
  bool IsTooFarFromBoundaries(TVector3 const& point) const;
  TVector3 PretendAtBoundary(TVector3 const& point) const;

  TVector3 ElectronDiverterPosOffsets(TVector3 const& point) const;
  bool   fEnableElectronDiverterDistortions1;
  bool   fEnableElectronDiverterDistortions2;
  double fEDZCenterMin;
  double fEDZCenterMax;
  double fEDAXPosOffs;
  double fEDBZPosOffs;
  double fEDs;
  double fEDChargeLossZLowMin;
  double fEDChargeLossZLowMax;
  double fEDChargeLossZHighMin;
  double fEDChargeLossZHighMax;

  TSpline3* MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const;
  double InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const;
  
  std::string fInputFilename;
  TFile* fInputFile;
    
  TH3F* hDx_cal_pos;
  TH3F* hDy_cal_pos;
  TH3F* hDz_cal_pos;
  TH3F* hEx_cal_pos;
  TH3F* hEy_cal_pos;
  TH3F* hEz_cal_pos;
  		   
  TH3F* hDx_cal_neg;
  TH3F* hDy_cal_neg;
  TH3F* hDz_cal_neg;
  TH3F* hEx_cal_neg;
  TH3F* hEy_cal_neg;
  TH3F* hEz_cal_neg;

  std::vector<std::vector<TSpline3*>> spline_dx_fwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dy_fwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dz_fwd_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dx_bkwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dy_bkwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dz_bkwd_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dEx_neg;
  std::vector<std::vector<TSpline3*>> spline_dEy_neg;
  std::vector<std::vector<TSpline3*>> spline_dEz_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dx_fwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dy_fwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dz_fwd_pos;
  
  std::vector<std::vector<TSpline3*>> spline_dx_bkwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dy_bkwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dz_bkwd_pos;
  
  std::vector<std::vector<TSpline3*>> spline_dEx_pos;
  std::vector<std::vector<TSpline3*>> spline_dEy_pos;
  std::vector<std::vector<TSpline3*>> spline_dEz_pos;

  //for systematic propagation: save a copy of the nominal histograms to easily reset SCE
  TH3F* nominal_hDx_cal_pos;
  TH3F* nominal_hDy_cal_pos;
  TH3F* nominal_hDz_cal_pos;
  TH3F* nominal_hEx_cal_pos;
  TH3F* nominal_hEy_cal_pos;
  TH3F* nominal_hEz_cal_pos;

  TH3F* nominal_hDx_cal_neg;
  TH3F* nominal_hDy_cal_neg;
  TH3F* nominal_hDz_cal_neg;
  TH3F* nominal_hEx_cal_neg;
  TH3F* nominal_hEy_cal_neg;
  TH3F* nominal_hEz_cal_neg;

  bool _IsVaried;
};
  
#endif
