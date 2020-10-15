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

  // From PionAnalizer_module otput
  Int_t trigger_offset(detinfo::DetectorClocksData const& data)
  {
    return 500;//data.TPCClock().Ticks(data.TriggerOffsetTPC() * -1.);
  }

  // From PionAnalizer_module otput
  double sampling_rate(detinfo::DetectorClocksData const& data)
  {
    return 500;//data.TPCClock().TickPeriod() * 1.e3;
  }
  
  //--------------------------------------------------------------------
  CalorimetryAlg::CalorimetryAlg()
    : //fCalAmpConstants{std::vector<double>()}
    //, fCalAreaConstants{std::vector<double>()}
     fUseModBox{false}
    , fLifeTimeForm{0}
    , fDoLifeTimeCorrection{false}
  {


    fCalAmpConstants.push_back(0.582554e-3);
    fCalAmpConstants.push_back(1.16594e-3);

    //    fCalAreaConstants.push_back(0.544391e-2);
    //    fCalAreaConstants.push_back(2.0376e-2);

    fCalAreaConstants.push_back(1e-3);
    fCalAreaConstants.push_back(1e-3);
    fCalAreaConstants.push_back(1e-3);

    std::cout << "anselmo -1: " << fCalAreaConstants[2]  << std::endl;
    
    if (fLifeTimeForm != 0 and fLifeTimeForm != 1) {
      /*
      throw cet::exception("CalorimetryAlg")
        << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
        << "Must select either '0' for exponential or '1' for exponential + "
           "constant.\n";
      */
    }
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           AnaHitPD const& hit,
                           double const pitch,
                           double const T0) const
  {
    return dEdx_AMP(
      clock_data, det_prop, hit.PeakAmplitude() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           double const dQ,
                           double const time,
                           double const pitch,
                           unsigned int const plane,
                           double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AMP(clock_data, det_prop, dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                           detinfo::DetectorPropertiesData const& det_prop,
                           double const dQdx,
                           double const time,
                           unsigned int const plane,
                           double const T0) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0);
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            AnaHitPD const& hit,
                            double const pitch,
                            double const T0) const
  {

    std::cout << "anselmo 3: " << hit.Integral()  << " " << pitch << " " <<  hit.PeakTime() << " " << hit.WireID().Plane << " " <<  T0 << std::endl;
    return dEdx_AREA(
      clock_data, det_prop, hit.Integral() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double const dQ,
                            double const time,
                            double const pitch,
                            unsigned int const plane,
                            double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AREA(clock_data, det_prop, dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double const dQdx,
                            double const time,
                            unsigned int const plane,
                            double const T0) const
  {

    double const fADCtoEl = fCalAreaConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm


    return dEdx_from_dQdx_e(clock_data, det_prop, dQdx_e, time, T0);
  }

  // Apply Lifetime and recombination correction.
  double
  CalorimetryAlg::dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                                   detinfo::DetectorPropertiesData const& det_prop,
                                   double dQdx_e,
                                   double const time,
                                   double const T0) const
  {
    if (fDoLifeTimeCorrection) {
      dQdx_e *= LifetimeCorrection(clock_data, det_prop, time, T0); // (dQdx_e in e/cm)
    }

    if (fUseModBox) { return ModBoxCorrection(dQdx_e); }

    return BirksCorrection(dQdx_e);
  }

  //------------------------------------------------------------------------------------//
  // for the time being copying from Calorimetry.cxx - should be decided where
  // to keep it.
  // ----------------------------------------------------------------------------------//
  double
  calo::CalorimetryAlg::LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
                                           detinfo::DetectorPropertiesData const& det_prop,
                                           double const time,
                                           double const T0) const
  {
    float const t = time - trigger_offset(clock_data);
    double const timetick = sampling_rate(clock_data) * 1.e-3; // time sample in microsec
    double const adjusted_time = t * timetick - T0 * 1e-3;     //  (in microsec)


    double ElectronLifetime = 35000; // from PionAnalyzer_module output
    
    assert(fLifeTimeForm < 2);
    if (fLifeTimeForm == 0) {
      // Exponential form
      double const tau = ElectronLifetime; //det_prop.ElectronLifetime();
      return exp(adjusted_time / tau);
    }
    /*
    // Exponential+constant form
    auto const& elifetime_provider =
      art::ServiceHandle<lariov::ElectronLifetimeService const>()->GetProvider();
    return elifetime_provider.Lifetime(adjusted_time);
    */
    return 0;
  }


double calo::CalorimetryAlg::ModBoxCorrection 	( 	double  	dQdx	) 	const
 {
   // Modified Box model correction has better behavior than the Birks
   // correction at high values of dQ/dx.


   double kModBoxA        = 0.930;    
   double kModBoxB        = 0.212;    
   double kGeVToElectrons = 4.237e7;
   double Density = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia)
   double Efield = pdAnaUtils::ComputeTotalEField(0,0,0);
   
   double const rho = Density;                          // LAr density in g/cm^3
   double Wion = 1000. / kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
   double const E_field = Efield; // Electric Field in the drift region in KV/cm
   double const Beta = kModBoxB / (rho * E_field);
   double Alpha = kModBoxA;
   double const dEdx = (exp(Beta * Wion * dQdx) - Alpha) / Beta;
   
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

   std::cout << "dedx = " << dEdx << std::endl;
   
   return dEdx;
 }
  
} // namespace
