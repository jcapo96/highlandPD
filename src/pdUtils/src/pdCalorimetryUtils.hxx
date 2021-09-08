//
// from CalorimetryAlg of  andrzej.szelc@yale.edu
//
#ifndef pdCalorimetryUtils_h
#define pdCalorimetryUtils_h

#include "pdDataClasses.hxx"
#include <vector>

class pdCalorimetryUtils {
public:
  pdCalorimetryUtils();
      
  Float_t ElectronsFromADCPeak(Float_t adc,  unsigned short plane)  const{return adc / fCalAmpConstants[plane];}  
  Float_t ElectronsFromADCArea(Float_t area, unsigned short plane) const{return area / fCalAreaConstants[plane];}

  Float_t LifetimeCorrection(Float_t time,Float_t T0 = 0) const;
  
  Float_t BirksCorrection  (Float_t dQdX) const;
  Float_t ModBoxCorrection (Float_t dQdX) const;    
  
  Float_t dEdx_from_dQdx_e(Float_t dQdx_e,Float_t time,Float_t T0 = 0) const;
  Float_t dEdx_from_dQdx_e(Float_t dQdx_e) const;
  
  Float_t dEdx_from_dQdx(Float_t dQdx,Float_t time,Float_t T0 = 0) const;
  Float_t dEdx_from_dQdx(Float_t dQdx) const;
  
  Float_t dQdx_e_from_dQdx(Float_t dQdx) const;
  
  void SetElectronLifetime(Float_t lifetime){fElectronLifetime=lifetime;}
  void SetModBoxA(Float_t ModBoxA){fModBoxA=ModBoxA;}
  void SetModBoxB(Float_t ModBoxB){fModBoxB=ModBoxB;}

  //---- OBSOLETE from CalorimetryAlg -----
  
  double dEdx_AMP(const AnaHitPD& hit,                   const double pitch,                     const double T0 = 0) const;
  double dEdx_AMP(const double dQ,     const double time,const double pitch, const UInt_t plane, const double T0 = 0) const;
  double dEdx_AMP(const double dQdx,   const double time,                    const UInt_t plane, const double T0 = 0) const;
  
  double dEdx_AREA(const AnaHitPD& hit,                  const double pitch,                     const double T0 = 0) const;
  double dEdx_AREA(const double dQ,    const double time,const double pitch, const UInt_t plane, const double T0 = 0) const;
  double dEdx_AREA(const double dQdx,  const double time,                    const UInt_t plane, const double T0 = 0) const;
  double dEdx_AREA(const double dQdx,  const double time,const TVector3& pos,const UInt_t plane, const double T0 = 0) const;
  //---------------------------------------
  
protected:
  
  std::vector<Float_t> fCalAmpConstants;
  std::vector<Float_t> fCalAreaConstants;
  bool fUseModBox;
  int fLifeTimeForm;
  bool fDoLifeTimeCorrection;
  
  Float_t fElectronLifetime;
  Float_t fModBoxA;
  Float_t fModBoxB;

  Float_t fSamplingRate;
  Float_t fTriggerOffset;

  bool debug;
}; 

#endif 
