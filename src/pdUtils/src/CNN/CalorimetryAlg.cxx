//  \file CalorimetryAlg.cxx
//
//  \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// andrzej.szelc@yale.edu
//

// LArSoft includes
#include "CalorimetryAlg.hxx"
#include <assert.h> 
/*
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeService.h"
*/

#include "pdAnalysisUtils.hxx" 




namespace calo {

bool debug=false;

  // From PionAnalizer_module otput
  Int_t trigger_offset()
  {
    return 500;//data.TPCClock().Ticks(data.TriggerOffsetTPC() * -1.);
  }

  // From PionAnalizer_module otput
  double sampling_rate()
  {
    return 500;//data.TPCClock().TickPeriod() * 1.e3;
  }
  
  //--------------------------------------------------------------------
  CalorimetryAlg::CalorimetryAlg()
    : //fCalAmpConstants{std::vector<double>()}
    //, fCalAreaConstants{std::vector<double>()}
     fUseModBox{true}
    , fLifeTimeForm{0}
    , fDoLifeTimeCorrection{true}
  {


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
      throw cet::exception("CalorimetryAlg")
        << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
        << "Must select either '0' for exponential or '1' for exponential + "
           "constant.\n";
      */
    }

    fElectronLifetime = 35000; // from PionAnalyzer_module output
    fModBoxA          = 0.930;    
    fModBoxB          = 0.212;             

  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(AnaHitPD const& hit,
                           double const pitch,
                           double const T0) const
  {
    return dEdx_AMP(hit.PeakAmplitude / pitch, hit.PeakTime, hit.WireID.Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(double const dQ,
                           double const time,
                           double const pitch,
                           unsigned int const plane,
                           double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AMP(dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(double const dQdx,
                           double const time,
                           unsigned int const plane,
                           double const T0) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(dQdx_e, time, T0);
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse
  // ----------------------------------------------------------------------------------//

  double
  CalorimetryAlg::dEdx_AREA(AnaHitPD const& hit,
                            double const pitch,
                            double const T0) const
  {

    if (debug)
      std::cout << "- dEdx_AREA (1): integral = " << hit.Integral  << ", pitch = " << pitch << ", peakTime = " <<  hit.PeakTime << ", plane = " << hit.WireID.Plane << ", T0 = " <<  T0 << std::endl;
    return dEdx_AREA(hit.Integral / pitch, hit.PeakTime, hit.Position, hit.WireID.Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(double const dQ,
                            double const time,
                            double const pitch,
                            unsigned int const plane,
                            double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AREA(dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(double const dQdx,
                            double const time,
                            const TVector3& pos,
                            unsigned int const plane,
                            double const T0) const
  {

    double const fADCtoEl = 0.00615;//fCalAreaConstants[plane];
    //    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm


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

  // Apply Lifetime and recombination correction.
  double
  CalorimetryAlg::dEdx_from_dQdx_e(double dQdx_e,
                                   double const time,
                                   double const T0) const
  {
    if (fDoLifeTimeCorrection) {
      dQdx_e *= LifetimeCorrection(time, T0); // (dQdx_e in e/cm)
    }

    double dEdx = dEdx_from_dQdx_e(dQdx_e);

    if (debug)
      std::cout << "    - dEdx_from_dQdx_e: dQdx_e = " << dQdx_e << " dQdx_e (lifetime corr)= " << dQdx_e << " dEdx = " << dEdx << std::endl;
    
    return dEdx;
  }

  // Apply Lifetime and recombination correction.
  double
  CalorimetryAlg::dEdx_from_dQdx(double dQdx,
                                 double const time,
                                 double const T0) const
  {
    
    double dQdx_e = dQdx_e_from_dQdx(dQdx);

    if (debug)
      std::cout << "    - dEdx_from_dQdx: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;
    
    return dEdx_from_dQdx_e(dQdx_e,time,T0);
  }


  // Apply recombination correction only
  double
  CalorimetryAlg::dEdx_from_dQdx(double dQdx) const
  {

    double dQdx_e = dQdx_e_from_dQdx(dQdx);
    double dEdx = dEdx_from_dQdx_e(dQdx_e);

    if (debug)
      std::cout << "    - dEdx_from_dQdx: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << " dEdx = " << dEdx << std::endl;

    return dEdx;
  }


  // Apply recombination correction only
  double
  CalorimetryAlg::dEdx_from_dQdx_e(double dQdx_e) const
  {
    //    if (debug)
    //      std::cout << "    - dEdx_from_dQdx_e: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;

    if (fUseModBox) { return ModBoxCorrection(dQdx_e); }
    return BirksCorrection(dQdx_e);
  }


  // Apply recombination correction only
  double
  CalorimetryAlg::dQdx_e_from_dQdx(double dQdx) const
  {

    //    double calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //

    // calib factor taken from
    // https://cdcvs.fnal.gov/redmine/projects/protoduneana/repository/revisions/develop/entry/protoduneana/Utilities/ProtoDUNECalibration.fcl
    // and explained here: https://wiki.dunescience.org/wiki/DQdx_and_dEdx_calibration_instructions#dE.2Fdx_calibration_factor
    double calib_factor[3] = {1.011e-3, 1.011e-3, 1.011e-3};
    double norm_factor[3] = {1.0078, 1.0082, 0.9946};

    double dQdx_e = dQdx/calib_factor[2]*norm_factor[2];
    
    if (debug)
      std::cout << "    - dQdx_from_dQdx_e: dQdx (ADC) =  " << dQdx << " dQdx (e) = " << dQdx_e << std::endl;

    return dQdx_e;
  }
  
  //------------------------------------------------------------------------------------//
  // for the time being copying from Calorimetry.cxx - should be decided where
  // to keep it.
  // ----------------------------------------------------------------------------------//
  double
  calo::CalorimetryAlg::LifetimeCorrection(double const time,
                                           double const T0) const
  {
    float const t = time - trigger_offset();
    double const timetick = sampling_rate() * 1.e-3; // time sample in microsec
    double const adjusted_time = t * timetick - T0 * 1e-3;     //  (in microsec)
    
    assert(fLifeTimeForm < 2);
    if (fLifeTimeForm == 0) {
      // Exponential form
      double const tau = fElectronLifetime; //det_prop.ElectronLifetime();
      double corr = exp(adjusted_time / tau);
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


double calo::CalorimetryAlg::ModBoxCorrection(double dQdx_e) 	const
 {


   // 

   /*  https://arxiv.org/pdf/1306.1712.pdf
       
       dQ/dx = R_c x W_ion x dE/dx 




       dE/dx= (exp(βWion·(dQ/dx))−α)/β.

    */
   
   // Modified Box model correction has better behavior than the Birks
   // correction at high values of dQ/dx.
   
   // dq/dx should be in e/cm
   // dE/dx is returned in MeV/cm
   
   
   double kModBoxA        = fModBoxA;    
   double kModBoxB        = fModBoxB;
   //   double kGeVToElectrons = 4.237e7/5.2;  //anselmo
   double kGeVToElectrons = 4.237e7;
   double Density = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia)
   double Efield = pdAnaUtils::ComputeTotalEField(0,0,0);
   
   double const rho = Density;                          // LAr density in g/cm^3
   double Wion = 1000. / kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
   double const E_field = Efield; // Electric Field in the drift region in KV/cm
   double const Beta = kModBoxB / (rho * E_field);
   double Alpha = kModBoxA;
   double const dEdx = (exp(Beta * Wion * dQdx_e) - Alpha) / Beta;

   if (debug){
     std::cout << "  --> A,B  = " << kModBoxA << " " << kModBoxB << std::endl;
     std::cout << "  --> dQdx_e, dEdx (modbox)= " << dQdx_e << " " << dEdx << std::endl;
   }
   
   return dEdx;
 }

double calo::CalorimetryAlg::BirksCorrection 	( 	double  	dQdx	) 	const  
 {
   // Correction for charge quenching using parameterization from
   // S.Amoruso et al., NIM A 523 (2004) 275

   double kRecombA        = 0.800;    
   double kRecombk        = 0.0486;   
   double kGeVToElectrons = 4.237e7; 
   double Density = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
   double Efield = pdAnaUtils::ComputeTotalEField(0,0,0);
   
   double A3t = kRecombA;
   double K3t = kRecombk;                           // in KV/cm*(g/cm^2)/MeV
   double const rho = Density;                          // LAr density in g/cm^3
   double Wion = 1000. / kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
   double const E_field = Efield; // Electric Field in the drift region in KV/cm
   K3t /= rho;                      // KV/MeV
   double const dEdx = dQdx / (A3t / Wion - K3t / E_field * dQdx); // MeV/cm

   if (debug)
     std::cout << "  --> dedx (birks) = " << dEdx << std::endl;
   
   return dEdx;
 }
  
} // namespace
