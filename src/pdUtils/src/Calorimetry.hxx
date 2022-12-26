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
  void Initialize(SpaceCharge* sce);
  void SetSCE(SpaceCharge* sce, bool remove = true);
  
  void ApplySCECorrection(AnaParticlePD* part) const;
  void ApplySCECorrection(AnaHitPD& hit) const;
  
  void ApplyLifetimeCorrection(AnaParticlePD* part) const;
  void ApplyLifetimeCorrection(AnaHitPD& hit) const;
  void UndoLifetimeCorrection(AnaParticlePD* part) const;
  void UndoLifetimeCorrection(AnaHitPD& hit) const;

  void ApplyNormalization(AnaHitPD &hit) const;
  void ApplyXCalibration(AnaHitPD &hit) const;
  void ApplyYZCalibration(AnaHitPD &hit) const;
  
  void ApplyRecombination(AnaParticlePD* part) const;
  void ApplyRecombination(AnaHitPD &hit) const;
    
  void CalibratedQdx(AnaParticlePD* part) const; //think about naming
  void CalibratedQdx(AnaHitPD &hit) const; //think about naming

  double GetNormalization(AnaHitPD &hit) const;
  double GetXCalibration(AnaHitPD &hit) const;
  double GetYZCalibration(AnaHitPD &hit) const;

  double GetLifetime() const {return _Lifetime;}
  void SetLifetime(double Lifetime){_Lifetime = Lifetime;}
  
  double GetModBoxA() const {return _ModBoxA;}
  double GetModBoxB() const {return _ModBoxB;}
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

  double _Lifetime;
  double _Vdrift;
  double _APAx;

  double _CalAreaConstants[3];

  SpaceCharge* _sce;

};
  
#endif
