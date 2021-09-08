#include "pdCalorimetryUtils.hxx"
#include "pdAnalysisUtils.hxx" 
#include <assert.h> 
  
//*********************************************
pdCalorimetryUtils::pdCalorimetryUtils(){
//*********************************************    

  debug=false;
  
  fUseModBox=true;
  fLifeTimeForm=0;
  fDoLifeTimeCorrection=true;
  
  fCalAmpConstants.push_back(0.582554e-3);
  fCalAmpConstants.push_back(1.16594e-3);
  fCalAmpConstants.push_back(1.16594e-3);
  
  //    fCalAreaConstants.push_back(0.544391e-2);
  //    fCalAreaConstants.push_back(2.0376e-2);
  
  fCalAreaConstants.push_back(1e-3);
  fCalAreaConstants.push_back(1e-3);
  fCalAreaConstants.push_back(1e-3);
  
  if (fLifeTimeForm != 0 and fLifeTimeForm != 1) {
    /*
      throw cet::exception("pdCalorimetryUtils")
      << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
      << "Must select either '0' for exponential or '1' for exponential + "
      "constant.\n";
    */
  }
  
  fElectronLifetime = 35000; // from PionAnalyzer_module output
  fModBoxA          = 0.930;    
  fModBoxB          = 0.212;             
  
}


//*********************************************
Float_t pdCalorimetryUtils::dEdx_from_dQdx_e(Float_t dQdx_e,const Float_t time,const Float_t T0) const{
//*********************************************

  // Apply Lifetime and recombination correction.
  
  if (fDoLifeTimeCorrection) {
    dQdx_e *= LifetimeCorrection(time, T0); // (dQdx_e in e/cm)
  }
  
  Float_t dEdx = dEdx_from_dQdx_e(dQdx_e);
  
  if (debug)
    std::cout << "    - dEdx_from_dQdx_e: dQdx_e = " << dQdx_e << " dQdx_e (lifetime corr)= " << dQdx_e << " dEdx = " << dEdx << std::endl;
  
  return dEdx;
}

//*********************************************
Float_t pdCalorimetryUtils::dEdx_from_dQdx(Float_t dQdx,const Float_t time,const Float_t T0) const{
//*********************************************

  // Apply recombination correction.
  
  Float_t dQdx_e = dQdx_e_from_dQdx(dQdx);
  
  if (debug)
    std::cout << "    - dEdx_from_dQdx: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;
  
    return dEdx_from_dQdx_e(dQdx_e,time,T0);
}

//*********************************************
Float_t pdCalorimetryUtils::dEdx_from_dQdx(Float_t dQdx) const{
//*********************************************
  
  // Apply recombination correction only
  
  Float_t dQdx_e = dQdx_e_from_dQdx(dQdx);
  Float_t dEdx = dEdx_from_dQdx_e(dQdx_e);
  
  if (debug)
    std::cout << "    - dEdx_from_dQdx: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << " dEdx = " << dEdx << std::endl;
  
  return dEdx;
}

//*********************************************
Float_t pdCalorimetryUtils::dEdx_from_dQdx_e(Float_t dQdx_e) const{
//*********************************************

  // Apply recombination correction only
  
  //    if (debug)
  //      std::cout << "    - dEdx_from_dQdx_e: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;
  
  if (fUseModBox) { return ModBoxCorrection(dQdx_e); }
  return BirksCorrection(dQdx_e);
}

//*********************************************
Float_t pdCalorimetryUtils::dQdx_e_from_dQdx(Float_t dQdx) const{
//*********************************************


  // Apply recombination correction only
  
  //    Float_t calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //
  
  // calib factor taken from
  // https://cdcvs.fnal.gov/redmine/projects/protoduneana/repository/revisions/develop/entry/protoduneana/Utilities/ProtoDUNECalibration.fcl
  // and explained here: https://wiki.dunescience.org/wiki/DQdx_and_dEdx_calibration_instructions#dE.2Fdx_calibration_factor
  Float_t calib_factor[3] = {1.011e-3, 1.011e-3, 1.011e-3};
  Float_t norm_factor[3] = {1.0078, 1.0082, 0.9946};
  
  Float_t dQdx_e = dQdx/calib_factor[2]*norm_factor[2];
  
  if (debug)
    std::cout << "    - dQdx_from_dQdx_e: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;
  
  return dQdx_e;
}

//------------------------------------------------------------------------------------//
// for the time being copying from Calorimetry.cxx - should be decided where
// to keep it.
// ----------------------------------------------------------------------------------//

//*********************************************
Float_t pdCalorimetryUtils::LifetimeCorrection(const Float_t time,const Float_t T0) const{
//*********************************************
  
  float const t = time - fTriggerOffset;
  const Float_t timetick = fSamplingRate * 1.e-3; // time sample in microsec
  const Float_t adjusted_time = t * timetick - T0 * 1e-3;     //  (in microsec)
  
  assert(fLifeTimeForm < 2);
  if (fLifeTimeForm == 0) {
    // Exponential form
    const Float_t tau = fElectronLifetime; //det_prop.ElectronLifetime();
    Float_t corr = exp(adjusted_time / tau);
    if (debug)
      std::cout << "    - lifetime corr = " << corr << std::endl;
    return corr;
  }
  /*
  // Exponential+constant form
  auto const& elifetime_provider =
  art::ServiceHandle<lariov::ElectronLifetimeService const>()->GetProvider();
  return elifetime_provider.Lifetime(adjusted_time);
  */
  return 0;
}

//*********************************************
Float_t pdCalorimetryUtils::ModBoxCorrection(Float_t dQdx_e) 	const{
//*********************************************

  /*  https://arxiv.org/pdf/1306.1712.pdf
      
      dQ/dx = R_c x W_ion x dE/dx                         
      dE/dx= (exp(βWion·(dQ/dx))−α)/β.      
  */
  
  // Modified Box model correction has better behavior than the Birks
  // correction at high values of dQ/dx.
  
  // dq/dx should be in e/cm
  // dE/dx is returned in MeV/cm
  
  
  Float_t kModBoxA        = fModBoxA;    
  Float_t kModBoxB        = fModBoxB;
  //   Float_t kGeVToElectrons = 4.237e7/5.2;  //anselmo
  Float_t kGeVToElectrons = 4.237e7;
  Float_t Density = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia)
  Float_t Efield = pdAnaUtils::ComputeTotalEField(0,0,0);
  
  const Float_t rho = Density;                          // LAr density in g/cm^3
  Float_t Wion = 1000. / kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
  const Float_t E_field = Efield; // Electric Field in the drift region in KV/cm
  const Float_t Beta = kModBoxB / (rho * E_field);
  Float_t Alpha = kModBoxA;
  const Float_t dEdx = (exp(Beta * Wion * dQdx_e) - Alpha) / Beta;
  
  if (debug){
    std::cout << "  --> A,B  = " << kModBoxA << " " << kModBoxB << std::endl;
    std::cout << "  --> dQdx_e, dEdx (modbox)= " << dQdx_e << " " << dEdx << std::endl;
  }
  
  return dEdx;
}

//*********************************************
Float_t pdCalorimetryUtils::BirksCorrection 	(Float_t dQdx)const  {
//*********************************************

  // Correction for charge quenching using parameterization from
  // S.Amoruso et al., NIM A 523 (2004) 275

  Float_t kRecombA        = 0.800;    
  Float_t kRecombk        = 0.0486;   
  Float_t kGeVToElectrons = 4.237e7; 
  Float_t Density = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
  Float_t Efield = pdAnaUtils::ComputeTotalEField(0,0,0);
  
  Float_t A3t = kRecombA;
  Float_t K3t = kRecombk;                           // in KV/cm*(g/cm^2)/MeV
  const Float_t rho = Density;                          // LAr density in g/cm^3
  Float_t Wion = 1000. / kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
  const Float_t E_field = Efield; // Electric Field in the drift region in KV/cm
  K3t /= rho;                      // KV/MeV
  const Float_t dEdx = dQdx / (A3t / Wion - K3t / E_field * dQdx); // MeV/cm
  
  if (debug)
    std::cout << "  --> dedx (birks) = " << dEdx << std::endl;
  
  return dEdx;
}

//------------------------------------------------------------------------------------//
// Functions to calculate the dEdX based on the AMPLITUDE of the pulse
// ----------------------------------------------------------------------------------//
//*********************************************
double pdCalorimetryUtils::dEdx_AMP(const AnaHitPD& hit,const double pitch,const double T0) const{
//*********************************************
  return dEdx_AMP(hit.PeakAmplitude / pitch, hit.PeakTime, hit.WireID.Plane, T0);
}

//*********************************************
double pdCalorimetryUtils::dEdx_AMP(const double dQ,const double time,const double pitch,const UInt_t plane,const double T0) const{
//*********************************************

  const double dQdx = dQ / pitch; // in ADC/cm
  return dEdx_AMP(dQdx, time, plane, T0);
}

//*********************************************
double pdCalorimetryUtils::dEdx_AMP(const double dQdx,const double time,const UInt_t plane,const double T0) const{
//*********************************************
  const double fADCtoEl = fCalAmpConstants[plane];
  const double dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
  return dEdx_from_dQdx_e(dQdx_e, time, T0);
}

//------------------------------------------------------------------------------------//
// Functions to calculate the dEdX based on the AREA of the pulse
// ----------------------------------------------------------------------------------//

//*********************************************
double pdCalorimetryUtils::dEdx_AREA(const AnaHitPD& hit,const double pitch,const double T0) const{
//*********************************************
  
  if (debug)
    std::cout << "- dEdx_AREA (1): integral = " << hit.Integral  << ", pitch = " << pitch << ", peakTime = " <<  hit.PeakTime << ", plane = " << hit.WireID.Plane << ", T0 = " <<  T0 << std::endl;
  return dEdx_AREA(hit.Integral / pitch, hit.PeakTime, hit.Position, hit.WireID.Plane, T0);
}

//*********************************************
double pdCalorimetryUtils::dEdx_AREA(const double dQ,const double time,const double pitch,const UInt_t plane,const double T0) const{
//*********************************************    

  const double dQdx = dQ / pitch; // in ADC/cm
  return dEdx_AREA(dQdx, time, plane, T0);
}

//*********************************************
double pdCalorimetryUtils::dEdx_AREA(const double dQdx,const double time,const UInt_t plane,const double T0) const{
//*********************************************

  TVector3 pos(0,0,0);
  return dEdx_AREA(dQdx, time, pos, plane, T0);
}

//*********************************************
double pdCalorimetryUtils::dEdx_AREA(const double dQdx,const double time,const TVector3& pos,const UInt_t plane,const double T0) const{
//*********************************************
  
  const double fADCtoEl = 0.00615;//fCalAreaConstants[plane];   //TODO
  //    const double dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
  
  //    double calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //
  //  double norm_factor[3] = {1.0078, 1.0082, 0.9947};
  //    double norm_factor[3] = {1.0078, 1.0082, 0.9946};
  
  
  // dq/dx should be in e/cm
  // dE/dx is returned in MeV/cm
  
  double dQdx_e = dQdx;///calib_factor[2]*norm_factor[2];
  
  
  Float_t cal_dQdx_e = pdAnaUtils::ComputeCalibrateddQdX(dQdx_e, pos);
  
  if (debug)
    std::cout << "  - dEdx_AREA (2): dQdx = " << dQdx  << ", fADCtoEl = " << fADCtoEl << ", dQdx_e = " <<  dQdx_e << ", cal_dQdx = " << cal_dQdx_e << std::endl;
  
  return dEdx_from_dQdx_e(cal_dQdx_e, time, T0);
}

  

