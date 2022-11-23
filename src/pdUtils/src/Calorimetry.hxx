//Class with calorimetry calibration

#ifndef Calorimetry_h
#define Calorimetry_h

#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TVector3.h"

#include "pdDataClasses.hxx"
#include "SpaceCharge.hxx"


class Calorimetry{
public:
  
  Calorimetry();
  virtual ~Calorimetry();

  void Initialize();

  void CalibratedQdx(AnaHitPD &hit) const; //think about naming

  void ApplyNormalization(AnaHitPD &hit) const;
  void ApplyXCalibration(AnaHitPD &hit) const;
  void ApplyYZCalibration(AnaHitPD &hit) const;
  void ApplyRecombination(AnaHitPD &hit) const;
  void ApplyRecombination(AnaParticlePD* part) const;
    
  double GetNormalization(AnaHitPD &hit) const;
  double GetXCalibration(AnaHitPD &hit) const;
  double GetYZCalibration(AnaHitPD &hit) const;

  double GetModBoxA(){return _ModBoxA;}
  double GetModBoxB(){return _ModBoxB;}
  void SetModBoxA(double ModBoxA){_ModBoxA = ModBoxA;}
  void SetModBoxB(double ModBoxB){_ModBoxB = ModBoxB;}
  void ResetModBoxParameters();

protected:
  
  void CreateNormHistogram();
  void CreateXCalHistogram();
  void CreateYZCalHistogram();

  double GetElectronsFromADCArea(double dQdx, int PlaneID) const {return dQdx/_CalAreaConstants[PlaneID];}

  TH1F* _h_norm;
  TH1F* _h_XCal[3];
  TH2F* _h_YZCal[3][2];

  double _Efield;
  double _ModBoxA;
  double _ModBoxB;
  
  double _CalAreaConstants[3];

  SpaceCharge* _sce;

};
  
#endif
