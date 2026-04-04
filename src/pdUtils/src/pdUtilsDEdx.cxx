#include "pdMomReconstruction.hxx"
#include "pdUtilsDEdx.hxx"
#include "pdUtilsTrack.hxx"
#include <TF1.h>
#include <TMath.h>
#include <Math/VavilovAccurate.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TFile.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern TProfile* PionTemplate;

//***************************************************************
double pdAnaUtils::GetDensityCorrection(double beta, double gamma){
//***************************************************************

  const double density_C  = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a  = 0.19559;
  const double density_k  = 3.0;

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

  const double rho = 1.39;
  const double K   = 0.307075;
  const double A   = 39.948;
  const double I   = 188.0e-6;
  const double Me  = 0.511;

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = GetDensityCorrection(beta, gamma);

  double f = rho * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0);

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

  const double Me  = 0.511;

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));

  return Wmax;
}

//***************************************************************
double pdAnaUtils::GetLandauXi(double KE, double dx, double mass){
//***************************************************************

  const double rho = 1.39;
  const double K   = 0.307075;
  const double A   = 39.948;

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

  if(par[0] < 0.01){
    this_vav = TMath::Landau(y);
    this_vav =this_vav / a;
  }
  else if(par[0] > 10.){
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =TMath::Gaus(y, mu, sigma);
  }
  else{
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
    out.momentumGeV = static_cast<Float_t>(pdMomShared::RangeCmToMomentumGeV(R_eff, 211, tg_ke, mass));
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
    double range = tg->GetPointX(i);
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

  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  std::string ssparticle;
  Float_t mass;
  if (!PdgToMassAndKeName(PDG, mass, ssparticle)) return -999.;
  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return -999.;

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);
  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), 2, dedx, rr, landauTruncMinRRCm,
                             landauTailHitDropFraction)) {
    return -999.;
  }
  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;

  return result;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                            double landauTruncMinRRCm, double landauTailHitDropFraction){
//***************************************************************

  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  std::string ssparticle;
  Float_t mass;
  if (!PdgToMassAndKeName(PDG, mass, ssparticle)) return -999.;
  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return -999.;

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);
  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, maxRR, 2, dedx, rr, landauTruncMinRRCm, landauTailHitDropFraction)) {
    return -999.;
  }
  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;

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
