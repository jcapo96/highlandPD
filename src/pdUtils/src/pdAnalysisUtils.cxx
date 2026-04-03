#include "pdAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "AnalysisUtils.hxx"
#include <TH3F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/VavilovAccurate.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TVectorD.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <set>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cstdio>
//data for range-momentum conversion, muons
//http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf divided by LAr density for cm

float Range_grampercm_new[29] = {
  9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396,
  1.063E2/1.396,  1.725E2/1.396, 2.385E2/1.396, 4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396,
  2.297E3/1.396,  4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396, 1.013E4/1.396, 1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396,
  4.326E4/1.396,  5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396};

float KE_MeV_new[29]= {
  10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
  400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
  20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

//data for range-momentum conversion, protons

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


TFile* dEdX_template_file = new TFile( (std::string(getenv("PDUTILSROOT"))+"/data/dEdxrestemplates.root").c_str(), "OPEN" );
std::map< int, TProfile* > templates;

TProfile* ProtonTemplate = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
TProfile* MuonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu" );
TProfile* KaonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_ka" );
TProfile* PionTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_pi" );
/*
templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
*/

//*****************************************************************************
Float_t pdAnaUtils::ComputeRangeMomentum(double trkrange, int pdg){
//*****************************************************************************

  /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
     website:
     http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf

     CSDA table values:
     float Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0,
     6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
     2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
     4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
     4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}; float KE_MeV[30] = {10, 14,
     20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000,
     4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000,
     200000, 300000, 400000};

     Functions below are obtained by fitting polynomial fits to KE_MeV vs
     Range (cm) graph. A better fit was obtained by splitting the graph into
     two: Below range<=200cm,a polynomial of power 4 was a good fit; above
     200cm, a polynomial of power 6 was a good fit

     Fit errors for future purposes:
     Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2
     err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7; Above 200cm,
     Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3
     err=5.9155E-9; p4 err=9.71709E-13; p5 err=7.22381E-17;p6
     err=1.9709E-21;*/

  //*********For muon, the calculations are valid up to 1.91E4 cm range
  //corresponding to a Muon KE of 40 GeV**********//

  /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
    website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

    CSDA values:
    double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
    1500, 2000, 2500, 3000, 4000, 5000};

    double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
    2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
    9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
    2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
    7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};

    Functions below are obtained by fitting power and polynomial fits to
    KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
    graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
    polynomial of power 6 was a good fit

    Fit errors for future purposes:
    For power function fit: a=0.388873; and b=0.00347075
    Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
    err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
    err=2.96958E-17;*/

  //*********For proton, the calculations are valid up to 3.022E3 cm range
  //corresponding to a Muon KE of 5 GeV**********//

  if (trkrange < 0 || std::isnan(trkrange)) {
    //    std::cout << "TrackMomentumCalculator   "
    //              << "Invalid track range " << trkrange << " return -1" << std::endl;
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
  } else
    KE = -999;

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
    //    std::cout << "CSDARangeCalculator   "
    //              << "Invalid beam mom " << beammom << " return -1" << std::endl;
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
  }
  else
    CSDARange=-1;

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

  // convert to GeV
  return kinetic/units::GeV;
}

//********************************************************************
void pdAnaUtils::ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance){
//********************************************************************

  distance.clear();
  if (!part) return;
  if (part->Daughters.empty()) return;

  for(size_t i = 0; i < part->Daughters.size(); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);

    //compute distance between part and daughter. Assume dau can be reconstructed backwards
    Float_t dis1 = sqrt(pow(part->PositionEnd[0]-dau->PositionStart[0],2)+pow(part->PositionEnd[1]-dau->PositionStart[1],2)+pow(part->PositionEnd[2]-dau->PositionStart[2],2));
    Float_t dis2 = sqrt(pow(part->PositionEnd[0]-dau->PositionEnd[0],2)+pow(part->PositionEnd[1]-dau->PositionEnd[1],2)+pow(part->PositionEnd[2]-dau->PositionEnd[2],2));

    //return lowest value
    if(dis1<dis2)distance.push_back(dis1);
    else distance.push_back(dis2);
  }
}

//********************************************************************
AnaTrueParticle* pdAnaUtils::GetBeamTrueParticle(const AnaSpillB& spill){
//********************************************************************

  AnaTrueParticle* beampart = NULL;

  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;

  int true_id = 0;
  if(beamPart)
    if(beamPart->TrueObject)
      true_id = static_cast<AnaTrueParticleB*>(beamPart->TrueObject)->ID;

  if(spill.TrueParticles.size() > 0){
    for(int i = 0; i < (int)spill.TrueParticles.size(); i++){
      if(true_id == spill.TrueParticles[i]->ID){
        beampart = static_cast<AnaTrueParticle*>(spill.TrueParticles[i]);
        break;
      }
    }
  }

  return beampart;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(AnaEventB* event, Int_t ID){
//*****************************************************************************

  if(!event) return NULL;

  // Get all reconstructed tracks in the event
  AnaTrueParticleB** trueParticles = event->TrueParticles;
  Int_t nTrueParts                 = event->nTrueParticles;

  if(!trueParticles || nTrueParts <= 0) return NULL;

  for (Int_t i=0;i<nTrueParts;i++){
    if(trueParticles[i] && trueParticles[i]->ID == ID){
	return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
      }
    }
  return NULL;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID){
//*****************************************************************************

  // Get all reconstructed tracks in the event
  for (UInt_t i=0;i<trueParticles.size();i++){
    if (trueParticles[i]->ID == ID){
      return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
    }
  }

  return NULL;
}

//*****************************************************************************
AnaParticlePD* pdAnaUtils::GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID){
//*****************************************************************************

  // loop over reconstructed tracks
  for(UInt_t i = 0; i < particles.size(); i++){
    AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(particles[i]->TrueObject);
    if(!truepart)continue;
    if(truepart->ID == true_ID)
      return static_cast<AnaParticlePD*>(particles[i]);
  }

  return NULL;
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID(const AnaParticlePD& part, const int pdg ){
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

  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
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


      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);

  return std::make_pair(pid_chi2, npt);
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR){
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

  //check particles' length is at least as long as maximum RR
  // if(part.Length<RR)
  //   return std::make_pair(9999., -1);
  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    //break whenever above RR upper limit, 26 cm is maximum
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


      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);

  return std::make_pair(pid_chi2, npt);
}

//*****************************************************************************
bool pdAnaUtils::isBeamLike(AnaParticlePD* part, AnaBeamPD* beam ){
//*****************************************************************************

  if (!part) return false;

  //From Owen Goodwins studies
  Float_t mccuts[7]  ={-3.,  7., -8.,  7., 27.5, 32.5, 0.93};
  Float_t datacuts[7]={ 0., 10., -5., 10., 30.,  35.0, 0.93};

  //cast the beam particle
  AnaParticlePD* beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if (!beampart) return false;

  Float_t beampos[3],beamdir[3], dist[3], dcos=0, cuts[7];

  // different way of obtaining the beam position and angle for DATA and MC
  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  // For MC
  if (trueBeamPart){
    for (int i=0;i<3;i++){
      beampos[i] = trueBeamPart->Position[i]-trueBeamPart->Position[2]*(trueBeamPart->Direction[i]/trueBeamPart->Direction[2]);
      beamdir[i] = trueBeamPart->Direction[i];
    }
    for (int i=0;i<7;i++) cuts[i] = mccuts[i];
  }
  else{
    // For Data
    if(beam->nMomenta != 1 || beam->nTracks != 1)return false;
    for (int i=0;i<3;i++){
      beampos[i] = beampart->PositionEnd[i];
      beamdir[i] = beampart->DirectionEnd[i];
    }
    for (int i=0;i<7;i++) cuts[i] = datacuts[i];
  }

  // compute the difference in position and cos(angle) considering that particle could have been reconstructed backwards
  if(part->PositionEnd[2]<part->PositionStart[2] && part->PositionEnd[2]!=-999){
    for (int i=0;i<3;i++){
      dist[i] = part->PositionEnd[i] - beampos[i];
      dcos   += -1*(part->DirectionEnd[i])*beamdir[i];
    }
  }
  else{
    for (int i=0;i<3;i++){
      dist[i] = part->PositionStart[i] - beampos[i];
      dcos   += part->DirectionStart[i]*beamdir[i];
    }
  }

  if(dist[0] < cuts[0] || dist[0] > cuts[1] || dist[1] < cuts[2] || dist[1] > cuts[3] || dist[2] < cuts[4] || dist[2] > cuts[5] || dcos < cuts[6]) return false;
  else return true;
}

//***************************************************************
AnaParticlePD* pdAnaUtils::GetBeamParticle(const AnaEventC& event){
//***************************************************************

  //Get the beam
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<const AnaEventB*>(&event)->Beam);
  if (!beam->BeamParticle) return NULL;

  // Get the beam particle
  return static_cast<AnaParticlePD*>(beam->BeamParticle);
}

//***************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueBeamParticle(const AnaEventC& event){
//***************************************************************

  // Get the beam particle
  AnaParticlePD* beampart = GetBeamParticle(event);
  if (!beampart) return NULL;

  // Get the true beam particle
  return static_cast<AnaTrueParticlePD*>(beampart->TrueObject);
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromHitPosition(const AnaParticlePD* part){
//***************************************************************

  Float_t length = 0.;

  //check if hit vector is empty
  if(part->Hits[2].empty()){
    //std::cout << "HitPosition vector empty! Returning -1" << std::endl;
    return -1;
  }

  // Initial position
  TVector3 disp(part->Hits[2][0].Position.X(),part->Hits[2][0].Position.Y(),part->Hits[2][0].Position.Z());

  // Add subsequent hits
  for(size_t i = 1; i < part->Hits[2].size(); ++i){
    if (part->Hits[2][i].Position.X() == -999) break;
    TVector3 pos(part->Hits[2][i].Position.X(),part->Hits[2][i].Position.Y(),part->Hits[2][i].Position.Z());
    disp -= pos;
    length += disp.Mag();
    disp = pos;
  }

  return length;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part){
//***************************************************************

  if(!part)return -999;

  Float_t length = 0.;

  int ntps = part->TrjPoints.size();
  int i0 = 0;
  double x0 = 0,y0 = 0,z0 = 0,dx = 0,dy = 0,dz = 0;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      x0 = part->TrjPoints[itp].Position.X();
      y0 = part->TrjPoints[itp].Position.Y();
      z0 = part->TrjPoints[itp].Position.Z();
      i0 = itp;
      break;
    }
  }

  for(int itp = i0; itp < ntps; itp++){
    if(!part->TrjPoints[itp].IsValid())continue;
    dx = part->TrjPoints[itp].Position.X()-x0;
    dy = part->TrjPoints[itp].Position.Y()-y0;
    dz = part->TrjPoints[itp].Position.Z()-z0;

    length += sqrt(dx*dx+dy*dy+dz*dz);

    x0 = part->TrjPoints[itp].Position.X();
    y0 = part->TrjPoints[itp].Position.Y();
    z0 = part->TrjPoints[itp].Position.Z();
  }

  return length;
}

//***************************************************************
void pdAnaUtils::ComputeParticlePositionAndDirection(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  int ntps = part->TrjPoints.size();
  if(ntps<1)return;

  int ifirst = -1;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      ifirst = itp;
      break;
    }
  }
  if(ifirst != -1){
    part->PositionStart[0] = part->TrjPoints[ifirst].Position.X();
    part->PositionStart[1] = part->TrjPoints[ifirst].Position.Y();
    part->PositionStart[2] = part->TrjPoints[ifirst].Position.Z();
    part->DirectionStart[0] = part->TrjPoints[ifirst].Direction.X();
    part->DirectionStart[1] = part->TrjPoints[ifirst].Direction.Y();
    part->DirectionStart[2] = part->TrjPoints[ifirst].Direction.Z();
  }

  int ilast  = -1;
  for(int itp = 1; itp < ntps; itp++){
    if(part->TrjPoints[ntps-itp].IsValid()){
      ilast = ntps-itp;
      break;
    }
  }
  if(ilast != -1){
    part->PositionEnd[0] = part->TrjPoints[ilast].Position.X();
    part->PositionEnd[1] = part->TrjPoints[ilast].Position.Y();
    part->PositionEnd[2] = part->TrjPoints[ilast].Position.Z();
    part->DirectionEnd[0] = part->TrjPoints[ilast].Direction.X();
    part->DirectionEnd[1] = part->TrjPoints[ilast].Direction.Y();
    part->DirectionEnd[2] = part->TrjPoints[ilast].Direction.Z();
  }
}

//***************************************************************
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> dEdx){

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(dEdx.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = dEdx.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + dEdx[i];
    counter ++;
  }

  return accumulated/counter;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD> hits){
//***************************************************************

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(hits.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = hits.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + hits.at(i).dEdx;
    counter ++;
  }

  return accumulated/counter;
}

//***************************************************************
Float_t pdAnaUtils::ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeDistanceMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(mother->PositionEnd[i] - daughter->PositionStart[i],2);
  dis = sqrt(dis);
  return dis;
}

//***************************************************************
Float_t pdAnaUtils::ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeCosMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + mother->DirectionEnd[i] * daughter->DirectionStart[i];
  return cos;
}

//***************************************************************
Float_t pdAnaUtils::ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange){
//***************************************************************

  if(!part){
    std::cout << "ComputeAveragedEdxOverResRange: particle does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  if(part->Hits[2].empty()){
    //std::cout << "ComputeAveragedEdxOverResRange: has no hits" << std::endl;
    //std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double sumdedx = 0;
  int nhits      = 0;
  for(size_t i = 0; i < part->Hits[2].size(); i++){
    // std::cout << "DEBUG: Residual range: " << part->Hits[2][i].ResidualRange << std::endl;
    // std::cout << "DEBUG: Dedx: " << part->Hits[2][i].dEdx << std::endl;
    // std::cout << "DEBUG: Dedx_calib: " << part->Hits[2][i].dEdx_calib << std::endl;
    // std::cout << "DEBUG: Dedx_SCE: " << part->Hits[2][i].dEdx_SCE << std::endl;
    // std::cout << "DEBUG: Dedx_NoSCE: " << part->Hits[2][i].dEdx_NoSCE << std::endl;
    if(part->Hits[2][i].ResidualRange < maxresrange){
      if(part->Hits[2][i].dEdx != -999 && part->Hits[2][i].dEdx < 1000){
        sumdedx += part->Hits[2][i].dEdx;
        nhits++;
      }
    }
  }

  return sumdedx/nhits;
}

//***************************************************************
bool pdAnaUtils::IsStoppingInFV(AnaParticlePD *part){
//***************************************************************

  if(!part)return false;

  bool ItIs = true;

  if((TMath::Abs(part->PositionStart[0])>350 ||
      part->PositionStart[1]<50 || part->PositionStart[1]>550 ||
      part->PositionStart[2]<50 || part->PositionStart[2]>645)
      &&
      (TMath::Abs(part->PositionEnd[0])>350 ||
       part->PositionEnd[1]<50 || part->PositionEnd[1]>550 ||
       part->PositionEnd[2]<50 || part->PositionEnd[2]>645))
    ItIs = false;

  return ItIs;
}

//***************************************************************
int pdAnaUtils::GetHitTPCid(AnaHitPD& hit){
//***************************************************************

  return GetPosTPCid(hit.Position);
}

//***************************************************************
int pdAnaUtils::GetPosTPCid(TVector3 pos){
//***************************************************************

  int TPCid = -1;

  if(pos.X() < 0){
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 1;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 5;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 9;
  }
  else{
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 2;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 6;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 10;
  }

  return TPCid;
}

//***************************************************************
void pdAnaUtils::EstimateHitsDirection(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  //loop over hits
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
    if(ihit < (int)part->Hits[2].size()-1){
      TVector3 dir = part->Hits[2][ihit+1].Position-part->Hits[2][ihit].Position;
      if(dir.Mag()!=0)dir.SetMag(1);
      part->Hits[2][ihit].Direction.SetXYZ(dir.X(),dir.Y(),dir.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(dir.X(),dir.Y(),dir.Z());
    }
    else{
      part->Hits[2][ihit].Direction.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
    }
  }
}

//***************************************************************
void pdAnaUtils::ComputeResidualRange(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  std::vector<double> delta; delta.clear();
  const int nh = (int)part->Hits[2].size();
  if (nh < 2) return;
  for (int ihit = 0; ihit < nh - 1; ihit++) {
    TVector3 diff = part->Hits[2][ihit + 1].Position - part->Hits[2][ihit].Position;
    delta.push_back(diff.Mag());
  }

  std::vector<double> new_rr;
  new_rr.reserve(nh);
  new_rr.push_back(delta[0] / 2.);
  for (int ihit = 1; ihit < nh; ++ihit) {
    new_rr.push_back(new_rr[ihit - 1] + delta[ihit - 1]);
  }

  for (int ihit = 0; ihit < nh; ihit++) part->Hits[2][ihit].ResidualRange = new_rr[ihit];
}

//***************************************************************
Double_t pdAnaUtils::ComputeDepositedEnergy(AnaParticlePD* part){
//***************************************************************

  Double_t E = -999;
  if(!part)return E;

  int nhits = part->Hits[2].size();
  if(nhits <= 0)return E;

  E = 0;
  for(int ihit = 1; ihit < nhits-1; ihit++){
    if(part->Hits[2][ihit].dEdx > 1000. || part->Hits[2][ihit].dEdx==-999 || part->Hits[2][ihit].Pitch < 0)
      continue;
    E += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }

  return E;
}

//***************************************************************
Double_t pdAnaUtils::EstimateTrueMomAtAPABorder(AnaParticlePD* part){
//***************************************************************

  Double_t momf = -999;
  if(!part)return momf;

  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return momf;

  if(part->PositionStart[2]>220)return momf;

  int nhits = part->Hits[2].size();
  if(nhits <= 0)return momf;

  double depE = 0;
  for(int ihit = 0; ihit < nhits; ihit++){
    if(part->Hits[2][ihit].Position.Z()>220)break;
    depE += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }

  double mass;
  if(abs(truePart->PDG) == 211)mass       = 139.57;
  else if(abs(truePart->PDG) == 321)mass  = 493.7;
  else if(abs(truePart->PDG) == 2212)mass = 938.27;
  else return momf;

  double momi = truePart->Momentum*1000;
  momf = sqrt(pow(sqrt(momi*momi+mass*mass)-depE,2)-mass*mass);

  return momf;
}

//*****************************************************************************
Float_t pdAnaUtils::ComputeCalorimetricMomentum(AnaParticlePD* part, int pdg, bool includeDecayProducts){
//*****************************************************************************

  if(!part) return -999.0;

  // Always use pion mass hypothesis for the main particle
  const double pionMass = 0.13957; // GeV/c^2

  int plane = 2; // Collection plane

  // Check if particle has hits
  if(part->Hits[plane].empty()) {
    // No hits available
    return -999.0;
  }

  // Calculate deposited kinetic energy from dE/dx integration
  // Note: framework base units are MeV for energy
  double kineticEnergy = 0.0;
  // std::cout << "DEBUG: Hits size: " << part->Hits[plane].size() << std::endl;
  for(size_t i = 0; i < part->Hits[plane].size(); i++){
    // Skip invalid hits
    // std::cout << "DEBUG: Dedx: " << part->Hits[plane][i].dEdx << std::endl;
    // std::cout << "DEBUG: Dedx_NoSCE: " << part->Hits[plane][i].dEdx_NoSCE << std::endl;
    // std::cout << "DEBUG: Dedx: " << part->Hits[plane][i].dEdx_SCE << std::endl;
    // std::cout << "DEBUG: Pitch: " << part->Hits[plane][i].Pitch << std::endl;
    if(part->Hits[plane][i].dEdx == -999 ||
       part->Hits[plane][i].Pitch < 0) {
      continue;
    }

    double dedx = part->Hits[plane][i].dEdx;
    double pitch = part->Hits[plane][i].Pitch;

    // std::cout << "DEBUG: Dedx: " << dedx << std::endl;
    // std::cout << "DEBUG: Pitch: " << pitch << std::endl;

    // Add energy deposited in this hit
    kineticEnergy += dedx * pitch;
  }

  // If requested, add energy from decay products (daughters)
  // Do NOT assume daughter mass - just add their deposited kinetic energy directly
  if(includeDecayProducts && !part->Daughters.empty()) {
    for(size_t i = 0; i < part->Daughters.size(); i++) {
      AnaParticlePD* daughter = static_cast<AnaParticlePD*>(part->Daughters[i]);
      if(!daughter) continue;

      // Calculate daughter's deposited kinetic energy without mass assumption
      if(!daughter->Hits[plane].empty()) {
        for(size_t j = 0; j < daughter->Hits[plane].size(); j++){
          // Skip invalid hits
          if(daughter->Hits[plane][j].dEdx == -999 ||
             daughter->Hits[plane][j].Pitch < 0) {
            continue;
          }

          double dedx = daughter->Hits[plane][j].dEdx;
          double pitch = daughter->Hits[plane][j].Pitch;

          // Add energy deposited in this daughter hit
          kineticEnergy += dedx * pitch;
        }
      }
    }
  }

  // Calculate total energy in MeV: E = KE + m
  // Note: pionMass is in GeV, convert to MeV
  double totalEnergy = kineticEnergy/1000 + pionMass; // GeV
  // std::cout << "DEBUG: Kinetic energy: " << kineticEnergy << std::endl;
  // std::cout << "DEBUG: Pion mass: " << pionMass << std::endl;
  // std::cout << "DEBUG: Total energy: " << totalEnergy << std::endl;

  // Calculate momentum in MeV: p = sqrt(E^2 - m^2)
  double momentumGeV = sqrt(totalEnergy*totalEnergy - pionMass*pionMass); // GeV
  // std::cout << "DEBUG: Momentum GeV: " << momentumGeV << std::endl;

  // Convert to GeV/c for consistency with framework
  double momentum = momentumGeV; // GeV

  return (Float_t)momentum;
}

//***************************************************************
Double_t pdAnaUtils::ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts){
//***************************************************************

  double distance = 9999;
  for(int ipart = 0; ipart < nparts; ipart++){
    AnaParticlePD* other = static_cast<AnaParticlePD*>(parts[ipart]);
    if(part->UniqueID == other->UniqueID)continue;
    double dis1 = 0;
    double dis2 = 0;
    for(int idis = 0; idis < 3; idis++){
      dis1 += pow(part->PositionEnd[idis]-other->PositionStart[idis],2);
      dis2 += pow(part->PositionEnd[idis]-other->PositionEnd[idis],2);
    }
    dis1 = sqrt(dis1);
    dis2 = sqrt(dis2);
    if(dis1 < distance)distance = dis1;
    if(dis2 < distance)distance = dis2;
  }

  return distance;
}


//***************************************************************
void pdAnaUtils::GetBeamQualityCuts(AnaEventPD* event,
				    double &mean_x, double &mean_y, double &mean_z,
				    double &sigma_x, double &sigma_y, double &sigma_z,
				    double &cos){
//***************************************************************

  //get nominal beam momentum. If none, set it to 1.
  AnaEventInfoPD* EventInfo = static_cast<AnaEventInfoPD*>(event->EventInfo);
  int NomBeamMom = (int)EventInfo->NominalBeamMom;
  if(NomBeamMom < 0 || NomBeamMom > 3)NomBeamMom = 1; //we still have no values for 6 and 7 GeV

  //get BQC parameters depending on beam mom and MC/data
  std::stringstream ssmom;
  ssmom << NomBeamMom;
  if(event->GetIsMC()){
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.MC."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
  else{
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.Data."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
}

//***************************************************************
double pdAnaUtils::GetDensityCorrection(double beta, double gamma){
//***************************************************************

  //Parameters for the density correction
  const double density_C  = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a  = 0.19559;
  const double density_k  = 3.0;

  //Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

//***************************************************************
double pdAnaUtils::GetdEdxBetheBloch(double KE, double mass){
//***************************************************************

  //Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar
  const double I   = 188.0e-6; // [MeV], mean excitation energy
  const double Me  = 0.511; // [Mev], mass of electron

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = GetDensityCorrection(beta, gamma);

  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}

//***************************************************************
double pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(TGraph* tg_ke, double range_cm) {
//***************************************************************
  constexpr double kKeFloorMeV = 1e-4;
  constexpr double kKeCapMeV = 1e9;
  if (!tg_ke || !std::isfinite(range_cm) || range_cm < 0.) return -1.;

  auto adopt_ke = [&](double k) -> double {
    if (!std::isfinite(k)) return -1.;
    if (k < kKeFloorMeV) k = kKeFloorMeV;
    if (k > kKeCapMeV) k = kKeCapMeV;
    return k;
  };

  if (range_cm > 0.) {
    const double ke_eval = tg_ke->Eval(range_cm);
    if (std::isfinite(ke_eval) && ke_eval >= kKeFloorMeV && ke_eval <= kKeCapMeV) return ke_eval;
  }

  std::vector<std::pair<double, double>> pts;
  pts.reserve(static_cast<size_t>(std::max(1, tg_ke->GetN())));
  for (int i = 0; i < tg_ke->GetN(); ++i) {
    double x = 0., y = 0.;
    tg_ke->GetPoint(i, x, y);
    if (x > 0. && std::isfinite(x) && std::isfinite(y)) pts.emplace_back(x, std::max(y, 0.));
  }
  if (pts.empty()) return -1.;
  std::sort(pts.begin(), pts.end());

  const double r = range_cm;
  if (r <= 0.) return kKeFloorMeV;

  if (r >= pts.back().first) {
    const double ke_hi = tg_ke->Eval(r);
    if (std::isfinite(ke_hi) && ke_hi >= kKeFloorMeV) return adopt_ke(ke_hi);
    return adopt_ke(pts.back().second);
  }

  if (r <= pts[0].first) {
    const double r0 = std::max(pts[0].first, 1e-12);
    const double k0 = pts[0].second;
    if (pts.size() >= 2) {
      const double r1 = std::max(pts[1].first, r0 + 1e-12);
      const double k1 = pts[1].second;
      const double slope = (k1 - k0) / (r1 - r0);
      double ke_ext = k0 + slope * (r - r0);
      if (!(ke_ext > 0.) || !std::isfinite(ke_ext)) ke_ext = k0 * (r / r0);
      return adopt_ke(ke_ext);
    }
    return adopt_ke(k0 * (r / r0));
  }

  for (size_t i = 1; i < pts.size(); ++i) {
    if (r <= pts[i].first + 1e-15) {
      const double r0 = pts[i - 1].first, k0 = pts[i - 1].second;
      const double r1 = pts[i].first, k1 = pts[i].second;
      if (std::abs(r1 - r0) < 1e-18) return adopt_ke(k0);
      const double t = (r - r0) / (r1 - r0);
      return adopt_ke(k0 + t * (k1 - k0));
    }
  }
  return adopt_ke(pts.back().second);
}

//***************************************************************
double pdAnaUtils::GetWmax(double KE, double mass){
//***************************************************************

  const double Me  = 0.511; // [Mev], mass of electron

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));

  return Wmax;
}

//***************************************************************
double pdAnaUtils::GetLandauXi(double KE, double dx, double mass){
//***************************************************************

  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  return xi;
}

//***************************************************************
double pdAnaUtils::dEdxPDF(double *x, double *par){
//***************************************************************

  ROOT::Math::VavilovAccurate vav;

  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav = vav.Pdf(y, par[0], par[1]);
    this_vav = this_vav / a;
  }

  return this_vav;
}

namespace {

bool PdgToMassAndKeName(Int_t PDG, Float_t& mass, std::string& ssparticle) {
  if (PDG == 13) {
    ssparticle = "muon";
    mass = 105.66f;
    return true;
  }
  if (PDG == 211) {
    ssparticle = "pion";
    mass = 139.57f;
    return true;
  }
  if (PDG == 321) {
    ssparticle = "kaon";
    mass = 493.677f;
    return true;
  }
  if (PDG == 2212) {
    ssparticle = "proton";
    mass = 938.272f;
    return true;
  }
  return false;
}

/// Among hits with RR >= truncMinRRCm, drop the largest-dEdx fraction `dropFrac` from the sample. Hits at smaller RR are never removed.
void DropLargestDedxFractionAmongHighRRHits(std::vector<double>& dedx, std::vector<double>& rr, double truncMinRRCm,
                                            double dropFrac) {
  if (dedx.size() != rr.size() || dedx.empty()) return;
  if (dropFrac <= 0. || truncMinRRCm <= 0.) return;

  std::vector<int> eligible;
  eligible.reserve(dedx.size());
  for (size_t i = 0; i < dedx.size(); ++i) {
    if (rr[i] >= truncMinRRCm) eligible.push_back(static_cast<int>(i));
  }
  if (eligible.size() < 2u) return;

  std::sort(eligible.begin(), eligible.end(), [&](int a, int b) {
    return dedx[static_cast<size_t>(a)] > dedx[static_cast<size_t>(b)];
  });
  const int nDrop = static_cast<int>(
      std::floor(dropFrac * static_cast<double>(eligible.size())));
  if (nDrop <= 0) return;

  std::unordered_set<int> dropIdx;
  for (int k = 0; k < nDrop && k < static_cast<int>(eligible.size()); ++k)
    dropIdx.insert(eligible[static_cast<size_t>(k)]);

  std::vector<double> nd, nr;
  nd.reserve(dedx.size());
  nr.reserve(rr.size());
  for (size_t i = 0; i < dedx.size(); ++i) {
    if (dropIdx.count(static_cast<int>(i))) continue;
    nd.push_back(dedx[i]);
    nr.push_back(rr[i]);
  }
  dedx.swap(nd);
  rr.swap(nr);
}

bool InteriorDedxRrSample(AnaParticlePD* part, double maxRR, int minPoints, std::vector<double>& dedx,
                          std::vector<double>& rr, double landauTruncMinRRCm, double landauTailHitDropFraction) {
  dedx.clear();
  rr.clear();
  if (!part || part->Hits[2].size() < 3u) return false;
  const int n = (int)part->Hits[2].size();
  const double kRRNoCap = 0.5 * std::numeric_limits<double>::max();
  const bool capRR = maxRR < kRRNoCap;
  for (int ihit = 1; ihit < n - 1; ++ihit) {
    const AnaHitPD& h = part->Hits[2][ihit];
    if (capRR && h.ResidualRange > maxRR) continue;
    dedx.push_back(static_cast<double>(h.dEdx));
    rr.push_back(static_cast<double>(h.ResidualRange));
  }
  DropLargestDedxFractionAmongHighRRHits(dedx, rr, landauTruncMinRRCm, landauTailHitDropFraction);
  return (int)dedx.size() >= minPoints;
}

double MeasuredTrackLengthCm(const AnaParticlePD* part, const std::vector<double>& rrInterior) {
  if (part && part->Length > 0.f && part->Length != -999.f) return static_cast<double>(part->Length);
  double mx = 0.;
  for (double r : rrInterior) {
    if (r > mx) mx = r;
  }
  return mx;
}

bool CollectionPlaneResidualRangeLooksUnset(const AnaParticlePD* part) {
  if (!part || part->Hits[2].empty()) return false;
  for (const AnaHitPD& h : part->Hits[2]) {
    if (h.ResidualRange > 0 && std::isfinite(static_cast<double>(h.ResidualRange))) return false;
  }
  return true;
}

// Load once: reopening ke_vs_range.root per fit was making annihilation vertex creation orders of magnitude slower.
TGraph* KeVsRangeGraphCached(const std::string& ssparticle) {
  static std::unordered_map<std::string, TGraph*> cache;
  static bool triedLoad = false;
  if (!triedLoad) {
    triedLoad = true;
    const char* root = getenv("PDUTILSROOT");
    if (root) {
      TFile* f = TFile::Open((std::string(root) + "/data/ke_vs_range.root").c_str(), "READ");
      if (f && !f->IsZombie()) {
        const char* names[] = {"muon", "pion", "kaon", "proton"};
        for (const char* nm : names) {
          auto* g = dynamic_cast<TGraph*>(f->Get(nm));
          if (!g) continue;
          auto* clone = dynamic_cast<TGraph*>(g->Clone());
          if (clone) cache[nm] = clone;
        }
        f->Close();
      }
      delete f;
    }
  }
  auto it = cache.find(ssparticle);
  return (it == cache.end()) ? nullptr : it->second;
}

pdAnaUtils::DEdxFreeRangeFitResult ParticleFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax, double step,
                                                       double maxRRForHits, int minInteriorPoints,
                                                       bool computeMomentum, double landauTruncMinRRCm,
                                                       double landauTailHitDropFraction) {
  pdAnaUtils::DEdxFreeRangeFitResult bad;
  if (!part || part->Hits[2].empty()) return bad;
  std::string ssparticle;
  Float_t mass = 0.f;
  if (!PdgToMassAndKeName(PDG, mass, ssparticle)) return bad;

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, maxRRForHits, minInteriorPoints, dedx, rr, landauTruncMinRRCm,
                             landauTailHitDropFraction))
    return bad;

  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return bad;

  TGraph* tg = new TGraph((int)dedx.size(), &rr[0], &dedx[0]);

  const double lenCm = MeasuredTrackLengthCm(part, rr);
  pdAnaUtils::DEdxFreeRangeFitResult result =
      pdAnaUtils::dEdxLikelihoodFreeRangeFit(tg, tg_ke, mass, 0., Lmax, step, lenCm, computeMomentum);

  delete tg;
  return result;
}

} // namespace

static TF1* FreeRangeLikelihoodPdf() {
  static TF1* pdf =
      new TF1("pdf_dedx_freerange_global", pdAnaUtils::dEdxPDF, -10., 20., 5);
  return pdf;
}

static pdAnaUtils::DEdxFreeRangeFitResult RunFreeRangeScan(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0, double Lmax,
                                                         double step, double measuredTrackLengthCm, bool computeMomentum) {
  pdAnaUtils::DEdxFreeRangeFitResult out;
  if (!tg || !tg_ke || tg->GetN() < 1 || step <= 0 || Lmax < L0) return out;

  const double width = 0.65;
  TF1* pdf = FreeRangeLikelihoodPdf();
  std::vector<double> L_v;
  std::vector<double> Likelihood_v;

  for (double L = L0; L <= Lmax + 1e-9; L += step) {
    double likelihood = 0.;
    for (int i = 0; i < tg->GetN(); i++) {
      const double range = tg->GetPointX(i) + L;
      const double dEdx = tg->GetPointY(i);
      if (!(range > 0.) || !std::isfinite(range)) continue;
      const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
      if (ke < 0. || !std::isfinite(ke)) continue;
      const double gamma = (ke / mass) + 1.0;
      const double beta = TMath::Sqrt(1. - (1.0 / (gamma * gamma)));
      const double xi = pdAnaUtils::GetLandauXi(ke, width, mass);
      const double Wmax = pdAnaUtils::GetWmax(ke, mass);
      const double kappa = xi / Wmax;
      const double dEdx_BB = pdAnaUtils::GetdEdxBetheBloch(ke, mass);
      const double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
      pdf->SetParameters(par);
      const double pval = pdf->Eval(dEdx);
      if (pval == 0.) continue;
      likelihood += std::log(pval);
    }
    L_v.push_back(L);
    Likelihood_v.push_back(likelihood);
  }

  if (Likelihood_v.empty()) return out;

  auto it = std::max_element(Likelihood_v.begin(), Likelihood_v.end());
  const int index = static_cast<int>(std::distance(Likelihood_v.begin(), it));
  out.logLikelihood = static_cast<Float_t>(Likelihood_v[index]);
  out.bestOffsetCm = static_cast<Float_t>(L_v[index]);

  if (computeMomentum && measuredTrackLengthCm > 0. && std::isfinite(measuredTrackLengthCm)) {
    const double R_eff = measuredTrackLengthCm + static_cast<double>(L_v[index]);
    const double KE = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, R_eff);
    if (KE > 0. && KE < 1e6 && std::isfinite(KE)) {
      const double pMeV = std::sqrt(KE * KE + 2.0 * static_cast<double>(mass) * KE);
      out.momentumGeV = static_cast<Float_t>(pMeV / 1000.0);
    }
  }
  return out;
}

//***************************************************************
Float_t pdAnaUtils::dEdxLikelihood(TGraph* tg, TGraph* tg_ke,
				   Float_t mass){
//***************************************************************

  double width = 0.65;
  TF1* pdf = new TF1("pdf", dEdxPDF, -10., 20., 5);
  double likelihood = 0;
  for(int i = 0; i < tg->GetN(); i++){
    double range = tg->GetPointX(i);//+L;
    double dEdx = tg->GetPointY(i);
    double ke = KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
    if (ke < 0. || !std::isfinite(ke)) continue;
    double gamma = (ke/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double xi = GetLandauXi(ke, width, mass);
    double Wmax = GetWmax(ke, mass);
    double kappa = xi / Wmax;
    double dEdx_BB = GetdEdxBetheBloch(ke, mass);
    double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
    pdf->SetParameters(par);
    if(pdf->Eval(dEdx) == 0)continue;
    else
      likelihood += log(pdf->Eval(dEdx));
  }
  delete pdf;
  return likelihood;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG, double landauTruncMinRRCm,
                                      double landauTailHitDropFraction){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);
  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), 2, dedx, rr, landauTruncMinRRCm,
                             landauTailHitDropFraction)) {
    file_ke->Close("R");
    return -999.;
  }
  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                            double landauTruncMinRRCm, double landauTailHitDropFraction){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);
  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, maxRR, 2, dedx, rr, landauTruncMinRRCm, landauTailHitDropFraction)) {
    file_ke->Close("R");
    return -999.;
  }
  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
pdAnaUtils::DEdxFreeRangeFitResult pdAnaUtils::dEdxLikelihoodFreeRangeFit(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0,
                                                            double Lmax, double step, double measuredTrackLengthCm,
                                                            bool computeMomentum){
//***************************************************************
  return RunFreeRangeScan(tg, tg_ke, mass, L0, Lmax, step, measuredTrackLengthCm, computeMomentum);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke,
					    Float_t mass){
//***************************************************************
  const pdAnaUtils::DEdxFreeRangeFitResult r =
      dEdxLikelihoodFreeRangeFit(tg, tg_ke, mass, 0., 10., 0.1, -1., false);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG,
                                                                  double landauTruncMinRRCm,
                                                                  double landauTailHitDropFraction){
//***************************************************************

  const pdAnaUtils::DEdxFreeRangeFitResult r =
      ParticleFreeRangeFit(part, PDG, 10., 0.1, std::numeric_limits<double>::max(), 1, false,
                           landauTruncMinRRCm, landauTailHitDropFraction);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
pdAnaUtils::DEdxFreeRangeFitResult pdAnaUtils::GetdEdxLikelihoodFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax,
                                                               double step, double landauTruncMinRRCm,
                                                               double landauTailHitDropFraction){
//***************************************************************
  return ParticleFreeRangeFit(part, PDG, Lmax, step, std::numeric_limits<double>::max(), 2, true,
                              landauTruncMinRRCm, landauTailHitDropFraction);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                                                        double landauTruncMinRRCm,
                                                                        double landauTailHitDropFraction){
//***************************************************************

  const pdAnaUtils::DEdxFreeRangeFitResult r =
      ParticleFreeRangeFit(part, PDG, 10., 0.1, maxRR, 1, false, landauTruncMinRRCm, landauTailHitDropFraction);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
pdAnaUtils::DEdxFreeRangeFitResult pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG,
                                                                      const double maxRR, double Lmax, double step,
                                                                      double landauTruncMinRRCm,
                                                                      double landauTailHitDropFraction){
//***************************************************************
  return ParticleFreeRangeFit(part, PDG, Lmax, step, maxRR, 2, true, landauTruncMinRRCm, landauTailHitDropFraction);
}

//***************************************************************
TMultiGraph* pdAnaUtils::MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax, double step,
                                                             double landauTruncMinRRCm, double landauTailHitDropFraction,
                                                             const char* xAxisTitle){
//***************************************************************
  if (!part || part->Hits[2].size() < 3u) return nullptr;
  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), 2, dedx, rr, landauTruncMinRRCm,
                             landauTailHitDropFraction))
    return nullptr;

  TGraph* tg_ke = KeVsRangeGraphCached("pion");
  if (!tg_ke) return nullptr;
  TGraph* tgFit = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());
  const double lenCm = MeasuredTrackLengthCm(part, rr);
  const DEdxFreeRangeFitResult fit = pdAnaUtils::dEdxLikelihoodFreeRangeFit(
      tgFit, tg_ke, 139.57f, 0., Lmax, step, lenCm, true);
  delete tgFit;
  if (!std::isfinite(static_cast<double>(fit.bestOffsetCm)) || fit.bestOffsetCm <= -998.f) return nullptr;

  const double L = static_cast<double>(fit.bestOffsetCm);
  std::vector<double> rrShift(static_cast<size_t>(dedx.size()));
  for (size_t i = 0; i < dedx.size(); ++i) {
    rrShift[i] = rr[i] + L;
  }

  auto* gMeas = new TGraph(static_cast<int>(rr.size()), rr.data(), dedx.data());
  gMeas->SetName("g_meas");
  gMeas->SetTitle("measured dE/dx vs RR");
  gMeas->SetMarkerColor(kBlack);
  gMeas->SetMarkerStyle(20);
  gMeas->SetMarkerSize(0.8f);

  auto* gShift = new TGraph(static_cast<int>(rrShift.size()), rrShift.data(), dedx.data());
  gShift->SetName("g_rr_shifted");
  gShift->SetTitle(Form("same dE/dx vs (RR+L_{best}), L_{best}=%.3f cm", L));
  gShift->SetMarkerColor(kRed);
  gShift->SetMarkerStyle(21);
  gShift->SetMarkerSize(0.8f);

  if (!PionTemplate) {
    delete gMeas;
    delete gShift;
    return nullptr;
  }

  std::vector<double> tplX, tplY;
  const int nb = PionTemplate->GetNbinsX();
  tplX.reserve(static_cast<size_t>(nb));
  tplY.reserve(static_cast<size_t>(nb));
  for (int b = 1; b <= nb; ++b) {
    if (PionTemplate->GetBinEntries(b) < 1) continue;
    const double x = PionTemplate->GetBinCenter(b);
    const double y = PionTemplate->GetBinContent(b);
    if (std::isfinite(x) && std::isfinite(y) && y > 0.) {
      tplX.push_back(x);
      tplY.push_back(y);
    }
  }
  if (tplX.empty()) {
    delete gMeas;
    delete gShift;
    return nullptr;
  }

  double maxTplR = 0.;
  double minTplR = std::numeric_limits<double>::infinity();
  for (double xv : tplX) {
    if (xv > maxTplR) maxTplR = xv;
    if (xv < minTplR) minTplR = xv;
  }

  constexpr double kRefCurveRrMaxCm = 500.;
  constexpr double kRefCurveExtendStepCm = 1.0;
  constexpr double kLowRStartCm = 0.02;
  constexpr double kPionMassMeV = 139.57;
  TGraph* tg_ke_pion = KeVsRangeGraphCached("pion");
  std::vector<std::pair<double, double>> tplPts;
  tplPts.reserve(tplX.size() + 550);
  for (size_t i = 0; i < tplX.size(); ++i) tplPts.emplace_back(tplX[i], tplY[i]);

  if (tg_ke_pion && std::isfinite(minTplR) && minTplR > kLowRStartCm + 1e-6) {
    for (double rr = kLowRStartCm; rr < minTplR - 0.5 * kRefCurveExtendStepCm + 1e-9;
         rr += kRefCurveExtendStepCm) {
      const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke_pion, rr);
      if (ke < 0. || !std::isfinite(ke)) continue;
      const double ybb = pdAnaUtils::GetdEdxBetheBloch(ke, static_cast<Float_t>(kPionMassMeV));
      if (!std::isfinite(ybb) || ybb <= 0.) continue;
      tplPts.emplace_back(rr, ybb);
    }
  }
  if (tg_ke_pion && maxTplR + kRefCurveExtendStepCm <= kRefCurveRrMaxCm + 1e-9) {
    for (double rr = maxTplR + kRefCurveExtendStepCm; rr <= kRefCurveRrMaxCm + 1e-9;
         rr += kRefCurveExtendStepCm) {
      const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke_pion, rr);
      if (ke < 0. || !std::isfinite(ke)) continue;
      const double ybb = pdAnaUtils::GetdEdxBetheBloch(ke, static_cast<Float_t>(kPionMassMeV));
      if (!std::isfinite(ybb) || ybb <= 0.) continue;
      tplPts.emplace_back(rr, ybb);
    }
  }
  std::sort(tplPts.begin(), tplPts.end());
  tplX.resize(tplPts.size());
  tplY.resize(tplPts.size());
  for (size_t i = 0; i < tplPts.size(); ++i) {
    tplX[i] = tplPts[i].first;
    tplY[i] = tplPts[i].second;
  }

  auto* gTpl = new TGraph(static_cast<int>(tplX.size()), tplX.data(), tplY.data());
  gTpl->SetName("g_pion_dedx_template");
  gTpl->SetTitle("pion ref: dedx_range_pi + dE/dx(BB) to 500 cm via KE(RR)");
  gTpl->SetLineColor(kBlue);
  gTpl->SetLineWidth(2);

  const char* xax = (xAxisTitle && xAxisTitle[0]) ? xAxisTitle : "Residual range [cm]";
  const std::string mgTitle = std::string(";") + xax + ";dE/dx [MeV/cm]";
  auto* mg = new TMultiGraph("mg_pion_freerange_dedx", mgTitle.c_str());
  mg->Add(gMeas, "P");
  mg->Add(gShift, "P");
  mg->Add(gTpl, "L");
  TH1F* frame = mg->GetHistogram();
  if (frame) {
    frame->SetTitle(mg->GetTitle());
    frame->GetXaxis()->SetTitle(xax);
    frame->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  }
  return mg;
}

//***************************************************************
// Usage example:
//   (legacy neutral/vertex-builder examples removed; use dedicated utils modules).
//***************************************************************
// Legacy duplicated neutral/vertex builders were removed.

//***************************************************************
// Usage example:
//   std::vector<double> fitParams;
//   pdAnaUtils::ExtrapolateTrack(particle, fitParams);
//   // fitParams: [x0, y0, z0, dx, dy, dz] for hits within 15 cm of DefinePosition
//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition,
                                  double trackFitDistanceFromStart){

  // Initialize output vector with 6 parameters
  fitParams.clear();
  fitParams.resize(6, -999.0);

  if (!part) {
    return;
  }

  // Get the reference position for distance calculations
  TVector3 referencePos = pdAnaUtils::DefinePosition(part, useStartPosition);
  if (referencePos.X() < -900) {
    return; // Invalid reference position
  }

  // Build reference travel direction.
  TVector3 travelDir(-999, -999, -999);
  if (useStartPosition) {
    travelDir.SetXYZ(part->DirectionStart[0], part->DirectionStart[1], part->DirectionStart[2]);
  } else {
    // If extrapolating from end, follow track backwards from end.
    travelDir.SetXYZ(-part->DirectionEnd[0], -part->DirectionEnd[1], -part->DirectionEnd[2]);
  }
  const bool hasValidTravelDir = (travelDir.Mag2() > 1e-10);
  if (hasValidTravelDir) travelDir = travelDir.Unit();

  // Build arc-length map from trajectory points (preferred for true along-track distance).
  std::vector<std::pair<TVector3, double>> trajectoryPointsWithDistance;
  if (part->TrjPoints.size() >= 2) {
    trajectoryPointsWithDistance.reserve(part->TrjPoints.size());
    double cumulative = 0.0;
    TVector3 prev;
    bool hasPrev = false;

    if (useStartPosition) {
      for (size_t i = 0; i < part->TrjPoints.size(); ++i) {
        const TVector3 pos = part->TrjPoints[i].Position;
        if (pos.Z() == -999) continue;
        if (hasPrev) cumulative += (pos - prev).Mag();
        trajectoryPointsWithDistance.push_back(std::make_pair(pos, cumulative));
        prev = pos;
        hasPrev = true;
      }
    } else {
      for (int i = static_cast<int>(part->TrjPoints.size()) - 1; i >= 0; --i) {
        const TVector3 pos = part->TrjPoints[i].Position;
        if (pos.Z() == -999) continue;
        if (hasPrev) cumulative += (pos - prev).Mag();
        trajectoryPointsWithDistance.push_back(std::make_pair(pos, cumulative));
        prev = pos;
        hasPrev = true;
      }
    }
  }

  auto computePathDistance = [&](const TVector3& position) -> double {
    if (!trajectoryPointsWithDistance.empty()) {
      double bestDist2 = 1e30;
      double bestArc = -1.0;
      for (const auto& tp : trajectoryPointsWithDistance) {
        const double d2 = (position - tp.first).Mag2();
        if (d2 < bestDist2) {
          bestDist2 = d2;
          bestArc = tp.second;
        }
      }
      return bestArc;
    }

    TVector3 delta = position - referencePos;
    return hasValidTravelDir ? delta.Dot(travelDir) : delta.Mag();
  };

  // Collect all 3D hit positions from all planes with their path distance from reference.
  std::vector<std::pair<TVector3, double>> hitPositionsWithDistance;

  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part->Hits[plane].size(); i++) {
      AnaHitPD& hit = part->Hits[plane][i];
      TVector3 position;

      // Use SCE-corrected position to match event display coordinate system
      if (hit.Position.Z() != -999) {
        position = hit.Position;
      } else {
        continue; // Skip invalid hits
      }

      const double pathDistance = computePathDistance(position);
      hitPositionsWithDistance.push_back(std::make_pair(position, pathDistance));
    }
  }

  // Need at least 2 points to fit a line
  if (hitPositionsWithDistance.size() < 2) {
    return;
  }

  const double fitWindowStart = std::max(0.0, trackFitDistanceFromStart);
  const double fitWindowEnd = fitWindowStart + std::max(0.0, trackLength);

  // Anchor the fitted line at the physical hit closest to fitWindowStart.
  TVector3 anchorPoint = referencePos;
  bool foundAnchor = false;
  double bestAnchorDelta = 1e30;
  for (const auto& hitPair : hitPositionsWithDistance) {
    if (hitPair.second < 0.0) continue;
    const double delta = std::abs(hitPair.second - fitWindowStart);
    if (!foundAnchor || delta < bestAnchorDelta) {
      bestAnchorDelta = delta;
      anchorPoint = hitPair.first;
      foundAnchor = true;
    }
  }
  if (!foundAnchor && hasValidTravelDir) {
    anchorPoint = referencePos + fitWindowStart * travelDir;
  }

  // Fit line to hits in the requested fit window along the travel direction.
  std::vector<TVector3> nearbyHits;
  for (const auto& hitPair : hitPositionsWithDistance) {
    if (hasValidTravelDir) {
      if (hitPair.second >= fitWindowStart && hitPair.second <= fitWindowEnd) {
        nearbyHits.push_back(hitPair.first);
      }
    } else if (hitPair.second >= fitWindowStart && hitPair.second <= fitWindowEnd) {
      nearbyHits.push_back(hitPair.first);
    }
  }

  if (nearbyHits.size() >= 2) {
    // Internal PCA line fit (kept inside ExtrapolateTrack to avoid exposing extra fit APIs).
    TVector3 centroid(0, 0, 0);
    for (const auto& pos : nearbyHits) {
      centroid += pos;
    }
    centroid *= (1.0 / nearbyHits.size());

    double covXX = 0, covXY = 0, covXZ = 0;
    double covYY = 0, covYZ = 0, covZZ = 0;
    for (const auto& pos : nearbyHits) {
      const double dx = pos.X() - centroid.X();
      const double dy = pos.Y() - centroid.Y();
      const double dz = pos.Z() - centroid.Z();
      covXX += dx * dx;
      covXY += dx * dy;
      covXZ += dx * dz;
      covYY += dy * dy;
      covYZ += dy * dz;
      covZZ += dz * dz;
    }

    const int nPoints = nearbyHits.size();
    covXX /= nPoints;
    covXY /= nPoints;
    covXZ /= nPoints;
    covYY /= nPoints;
    covYZ /= nPoints;
    covZZ /= nPoints;

    TMatrixD covMatrix(3, 3);
    covMatrix(0, 0) = covXX; covMatrix(0, 1) = covXY; covMatrix(0, 2) = covXZ;
    covMatrix(1, 0) = covXY; covMatrix(1, 1) = covYY; covMatrix(1, 2) = covYZ;
    covMatrix(2, 0) = covXZ; covMatrix(2, 1) = covYZ; covMatrix(2, 2) = covZZ;

    TMatrixDEigen eigen(covMatrix);
    TVectorD eigenValues = eigen.GetEigenValuesRe();
    TMatrixD eigenVectors = eigen.GetEigenVectors();

    int maxEigenIndex = 0;
    double maxEigenValue = eigenValues[0];
    for (int i = 1; i < 3; i++) {
      if (eigenValues[i] > maxEigenValue) {
        maxEigenValue = eigenValues[i];
        maxEigenIndex = i;
      }
    }

    TVector3 direction(eigenVectors(0, maxEigenIndex),
                       eigenVectors(1, maxEigenIndex),
                       eigenVectors(2, maxEigenIndex));
    direction = direction.Unit();

    fitParams[0] = anchorPoint.X();
    fitParams[1] = anchorPoint.Y();
    fitParams[2] = anchorPoint.Z();
    fitParams[3] = direction.X();
    fitParams[4] = direction.Y();
    fitParams[5] = direction.Z();
    // Enforce fitted direction to follow particle travel direction.
    if (hasValidTravelDir) {
      TVector3 fitDir(fitParams[3], fitParams[4], fitParams[5]);
      if (fitDir.Mag2() > 1e-10 && fitDir.Dot(travelDir) < 0.0) {
        fitParams[3] *= -1.0;
        fitParams[4] *= -1.0;
        fitParams[5] *= -1.0;
      }
    }
  }

}

//***************************************************************
// Define the position to use for calculations (distance, line fitting, etc.)
//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle, bool useStartPosition) {
//***************************************************************

  if (!particle) {
    return TVector3(-999, -999, -999);
  }

  // Choose between start and end position based on parameter
  if (useStartPosition) {
    // Use start position (for vertex particles)
  return TVector3(particle->PositionStart[0],
                  particle->PositionStart[1],
                  particle->PositionStart[2]);
  } else {
    // Use end position (for parent particles)
    return TVector3(particle->PositionEnd[0],
                    particle->PositionEnd[1],
                    particle->PositionEnd[2]);
  }

}


//***************************************************************
// Overloaded version for backward compatibility (defaults to start position)
//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle) {
//***************************************************************
  return DefinePosition(particle, true); // Default to start position for backward compatibility
}

//***************************************************************
// Overloaded version of ExtrapolateTrack for backward compatibility (defaults to start position)
//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength) {
//***************************************************************
  ExtrapolateTrack(part, fitParams, trackLength, true, 0.0); // Default to start position for backward compatibility
}

//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition){

  // Initialize output vector with 6 parameters
  fitParams.clear();
  fitParams.resize(6, -999.0);

  if (!part) {
    return;
  }

  // Get position and direction based on useStartPosition flag
  double x0, y0, z0, dx, dy, dz;

  if (useStartPosition) {
    x0 = part->Position[0];
    y0 = part->Position[1];
    z0 = part->Position[2];
    dx = part->Direction[0];
    dy = part->Direction[1];
    dz = part->Direction[2];
  } else {
    x0 = part->PositionEnd[0];
    y0 = part->PositionEnd[1];
    z0 = part->PositionEnd[2];
    dx = part->DirectionEnd[0];
    dy = part->DirectionEnd[1];
    dz = part->DirectionEnd[2];
  }

  // Check for valid values
  if (x0 == -999.0 || y0 == -999.0 || z0 == -999.0 ||
      dx == -999.0 || dy == -999.0 || dz == -999.0) {
    return;
  }

  // Normalize direction vector
  double norm = sqrt(dx*dx + dy*dy + dz*dz);
  if (norm > 0) {
    dx /= norm;
    dy /= norm;
    dz /= norm;
  } else {
    return; // Invalid direction
  }

  // Store the line parameters: [x0, y0, z0, dx, dy, dz]
  fitParams[0] = x0;
  fitParams[1] = y0;
  fitParams[2] = z0;
  fitParams[3] = dx;
  fitParams[4] = dy;
  fitParams[5] = dz;

}

//***************************************************************
// Overloaded version of ExtrapolateTrack for backward compatibility (defaults to start position)
//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaTrueParticlePD* part, std::vector<double>& fitParams, double trackLength) {
//***************************************************************
  ExtrapolateTrack(part, fitParams, trackLength, true); // Default to start position for backward compatibility
}

//***************************************************************
// Find the closest points between two 3D lines
//***************************************************************
double pdAnaUtils::FindClosestPointsBetweenLines(const std::vector<double>& line1Params,
                                                const std::vector<double>& line2Params,
                                                TVector3& closestPoint1,
                                                TVector3& closestPoint2) {

  // Check if both lines have valid parameters
  if (line1Params.size() != 6 || line2Params.size() != 6 ||
      line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    closestPoint1.SetXYZ(-999, -999, -999);
    closestPoint2.SetXYZ(-999, -999, -999);
    return -999.0;
  }

  // Extract line parameters
  // Line 1: P1(t1) = (x1, y1, z1) + t1 * (dx1, dy1, dz1)
  TVector3 point1(line1Params[0], line1Params[1], line1Params[2]);
  TVector3 dir1(line1Params[3], line1Params[4], line1Params[5]);

  // Line 2: P2(t2) = (x2, y2, z2) + t2 * (dx2, dy2, dz2)
  TVector3 point2(line2Params[0], line2Params[1], line2Params[2]);
  TVector3 dir2(line2Params[3], line2Params[4], line2Params[5]);

  // Normalize direction vectors
  dir1 = dir1.Unit();
  dir2 = dir2.Unit();

  // Vector between the two points
  TVector3 w0 = point1 - point2;

  // Calculate dot products
  double a = dir1.Dot(dir1);  // Should be 1 since dir1 is normalized
  double b = dir1.Dot(dir2);
  double c = dir2.Dot(dir2);  // Should be 1 since dir2 is normalized
  double d = dir1.Dot(w0);
  double e = dir2.Dot(w0);

  // Calculate denominator
  double denom = a * c - b * b;

  // Check if lines are parallel (denominator close to zero)
  if (fabs(denom) < 1e-10) {
    // Lines are parallel, use arbitrary point on line1 and project to line2
    closestPoint1 = point1;
    double t2 = e / c;  // Project point1 onto line2
    closestPoint2 = point2 + t2 * dir2;
  } else {
    // Lines are not parallel, find closest points
    double t1 = (b * e - c * d) / denom;
    double t2 = (a * e - b * d) / denom;

    closestPoint1 = point1 + t1 * dir1;
    closestPoint2 = point2 + t2 * dir2;
  }

  // Calculate the minimum distance
  double minDistance = (closestPoint1 - closestPoint2).Mag();

  return minDistance;
}

//********************************************************************
double pdAnaUtils::CalculateImpactParameter(const std::vector<double>& lineParams, const TVector3& point){
//********************************************************************

  if (lineParams.size() != 6 || lineParams[0] == -999.0) {
    return -999.0;
  }

  TVector3 linePoint(lineParams[0], lineParams[1], lineParams[2]);
  TVector3 lineDirection(lineParams[3], lineParams[4], lineParams[5]);
  lineDirection = lineDirection.Unit(); // Ensure it's normalized

  // Calculate distance from point to line
  TVector3 pointToLine = point - linePoint;
  TVector3 projection = (pointToLine.Dot(lineDirection)) * lineDirection;
  TVector3 perpendicular = pointToLine - projection;

  return perpendicular.Mag();
}

//*****************************************************************************
Float_t pdAnaUtils::ComputeTrueInvariantMass(const AnaTrueParticlePD& part1, const AnaTrueParticlePD& part2, Float_t mass1, Float_t mass2) {

  if (mass1 < 0. || mass2 < 0.)
    return -999.;

  // Calculate energies using E² = p² + m²
  Float_t E1 = sqrt(mass1*mass1 + part1.Momentum * part1.Momentum);
  Float_t E2 = sqrt(mass2*mass2 + part2.Momentum * part2.Momentum);

  // Calculate 3-momentum vectors from momentum magnitude and direction
  Float_t p1x = part1.Momentum * part1.Direction[0];
  Float_t p1y = part1.Momentum * part1.Direction[1];
  Float_t p1z = part1.Momentum * part1.Direction[2];

  Float_t p2x = part2.Momentum * part2.Direction[0];
  Float_t p2y = part2.Momentum * part2.Direction[1];
  Float_t p2z = part2.Momentum * part2.Direction[2];

  // Calculate invariant mass: M² = (E1 + E2)² - (p1 + p2)²
  Float_t totalE = E1 + E2;
  Float_t totalPx = p1x + p2x;
  Float_t totalPy = p1y + p2y;
  Float_t totalPz = p1z + p2z;
  Float_t totalP2 = totalPx * totalPx + totalPy * totalPy + totalPz * totalPz;

  Float_t M2 = totalE * totalE - totalP2;
  if (M2 > 0) {
    return sqrt(M2);
  }

  return -999.;
}