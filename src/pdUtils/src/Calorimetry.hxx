//Class with calorimetry calibration

#ifndef Calorimetry_h
#define Calorimetry_h

#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TVector3.h"

#include "pdDataClasses.hxx"


class Calorimetry{
public:
  
  Calorimetry();
  virtual ~Calorimetry();

  void Initialize();

  void CalibrateHit(AnaHitPD &hit) const; //think about naming
    
protected:
  
  void CreateNormHistogram();
  void CreateXcorrHistogram();
  void CreateYZcorrHistogram();

  void ApplyNormCorrection(AnaHitPD &hit) const;
  void ApplyXCorrection(AnaHitPD &hit) const;
  void ApplyYZCorrection(AnaHitPD &hit) const;
  void ApplyRecombination(AnaHitPD &hit) const;

  double GetElectronsFromADCArea(double dQdx, int PlaneID) const {return dQdx/_CalAreaConstants[PlaneID];}

  TH1F* _h_norm;
  TH1F* _h_Xcorr[3];
  TH2F* _h_YZcorr[3][2];

  double _Efield;
  double _ModBoxA;
  double _ModBoxB;
  
  double _CalAreaConstants[3];

};
  
#endif
