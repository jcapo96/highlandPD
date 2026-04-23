#include "pdUtilsRangePID.hxx"

#include "AnalysisUtils.hxx"
#include "TSpline.h"
#include <TFile.h>
#include <TGraph.h>
#include <TProfile.h>

#include <cmath>
#include <iostream>
#include <string>

// data for range-momentum conversion, muons
float Range_grampercm_new[29] = {
  9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396,
  1.063E2/1.396,  1.725E2/1.396, 2.385E2/1.396, 4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396,
  2.297E3/1.396,  4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396, 1.013E4/1.396, 1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396,
  4.326E4/1.396,  5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396};

float KE_MeV_new[29]= {
  10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
  400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
  20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

// data for range-momentum conversion, protons
double Range_gpercm_P_Nist[31]={
  1.887E-1/1.396, 3.823E-1/1.396, 6.335E-1/1.396, 1.296/1.396,   2.159/1.396,   7.375/1.396,   1.092E1/1.396, 2.215E1/1.396,
  3.627E1/1.396,  5.282E1/1.396,  7.144E1/1.396,  9.184E1/1.396, 1.138E2/1.396, 1.370E2/1.396, 1.614E2/1.396, 1.869E2/1.396,
  2.132E2/1.396,  2.403E2/1.396,  2.681E2/1.396,  2.965E2/1.396, 3.254E2/1.396, 3.548E2/1.396, 3.846E2/1.396, 4.148E2/1.396,
  4.454E2/1.396,  7.626E2/1.396,  1.090E3/1.396,  1.418E3/1.396, 1.745E3/1.396, 2.391E3/1.396, 3.022E3/1.396};

double KE_MeV_P_Nist[31]={
  10,   15,   20,   30,   40,   80,   100, 150,
  200,  250,  300,  350,  400,  450,  500, 550,
  600,  650,  700,  750,  800,  850,  900, 950,
  1000, 1500, 2000, 2500, 3000, 4000, 5000};

TGraph const KEvsR(29, Range_grampercm_new, KE_MeV_new);
TSpline3 const KEvsR_spline3("KEvsRS", &KEvsR);

TGraph const RvsKE(29, KE_MeV_new, Range_grampercm_new);
TSpline3 const RvsKE_spline3("RvsKES", &RvsKE);

TGraph const RvsKE_P(31, KE_MeV_P_Nist, Range_gpercm_P_Nist);
TSpline3 const RvsKE_P_spline3("RvsKE_P_S", &RvsKE_P);

TFile* dEdX_template_file = new TFile((std::string(getenv("PDUTILSROOT")) + "/data/dEdxrestemplates.root").c_str(), "OPEN");

TProfile* ProtonTemplate = (TProfile*)dEdX_template_file->Get("dedx_range_pro");
TProfile* MuonTemplate   = (TProfile*)dEdX_template_file->Get("dedx_range_mu");
TProfile* KaonTemplate   = (TProfile*)dEdX_template_file->Get("dedx_range_ka");
TProfile* PionTemplate   = (TProfile*)dEdX_template_file->Get("dedx_range_pi");

//*****************************************************************************
Float_t pdAnaUtils::ComputeRangeMomentum(double trkrange, int pdg){
//*****************************************************************************

  if (trkrange < 0 || std::isnan(trkrange)) {
    return -1.;
  }

  double KE, Momentum, M;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;

  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = KEvsR_spline3.Eval(trkrange);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    if (trkrange > 0 && trkrange <= 80)
      KE = 29.9317 * std::pow(trkrange, 0.586304);
    else if (trkrange > 80 && trkrange <= 3.022E3)
      KE =
        149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
        (4.34587E-6 * trkrange * trkrange * trkrange) +
        (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
        (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
        (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *
         trkrange);
    else
      KE = -999;
  } else {
    KE = -999;
  }

  if (KE < 0)
    Momentum = -999;
  else
    Momentum = std::sqrt((KE * KE) + (2 * M * KE));

  Momentum = Momentum / 1000;

  return Momentum;
}

//*****************************************************************************
Float_t pdAnaUtils::ComputeCSDARange(double beammom, int pdg){
//*****************************************************************************

  if (beammom < 0 || std::isnan(beammom)) {
    return -1.;
  }

  double KE, M, CSDARange;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;

  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_spline3.Eval(KE);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_P_spline3.Eval(KE);
  } else {
    CSDARange = -1;
  }

  return CSDARange;
}

//********************************************************************
Float_t pdAnaUtils::ComputeKineticEnergy(const AnaParticlePD &part) {
//********************************************************************

  int plane = 2;

  if(part.Hits[plane].size() < 1)return -1;

  double kinetic = 0;
  double res     = 0;

  for(size_t i = 0; i < part.Hits[plane].size(); i++){
    if(part.Hits[plane][i].dEdx_calib > 1000. || part.Hits[plane][i].dEdx_calib==-999)continue;
    double dedxi = part.Hits[plane][i].dEdx_calib;
    double Residualrangei = part.Hits[plane][i].ResidualRange;
    kinetic = kinetic + dedxi * fabs(Residualrangei - res);
    res = Residualrangei;
  }

  return kinetic/units::GeV;
}

//*****************************************************************************
std::pair<double, int> pdAnaUtils::Chi2PID(const AnaParticlePD& part, const int pdg ){
//*****************************************************************************

  double pid_chi2 = 0.;
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;

  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else if(pdg == 211)profile = PionTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }

  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);

  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    int bin = profile->FindBin( part.Hits[plane][i].ResidualRange );

    if( bin >= 1 && bin <= profile->GetNbinsX() ){

      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;
      }


      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;
      dedx_res *= part.Hits[plane][i].dEdx;

      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);

  return std::make_pair(pid_chi2, npt);
}

//*****************************************************************************
Float_t pdAnaUtils::Chi2PIDChi2PerHit(const AnaParticlePD* part, const int pdg){
//*****************************************************************************

  if(!part) return -999.f;
  const std::pair<double, int> r = Chi2PID(*part, pdg);
  if(r.first < 0.0 || r.second <= 0) return -999.f;
  return static_cast<Float_t>(r.first / static_cast<double>(r.second));
}

//*****************************************************************************
std::pair<double, int> pdAnaUtils::Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR){
//*****************************************************************************

  double pid_chi2 = 0.;
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;

  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }

  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);

  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    if(part.Hits[plane][i].ResidualRange > RR)continue;

    int bin = profile->FindBin( part.Hits[plane][i].ResidualRange );

    if( bin >= 1 && bin <= profile->GetNbinsX() ){

      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;
      }

      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;
      dedx_res *= part.Hits[plane][i].dEdx;

      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);

  return std::make_pair(pid_chi2, npt);
}
