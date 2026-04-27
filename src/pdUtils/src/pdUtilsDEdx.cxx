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

static TF1* FreeRangeLikelihoodPdf();
static bool BuildDedxPdfParams(double ke, Float_t mass, double width, double* par);

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

bool InteriorDedxRrSample(AnaParticlePD* part, double maxRR, int minPoints, std::vector<double>& dedx,
                          std::vector<double>& rr, int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                          double dedxMaxMeVcm) {
  dedx.clear();
  rr.clear();
  if (!part || part->Hits[2].empty()) return false;
  if (skipHitsFirst < 0) skipHitsFirst = 0;
  if (skipHitsLast < 0) skipHitsLast = 0;
  const int n = static_cast<int>(part->Hits[2].size());
  if (n < skipHitsFirst + skipHitsLast + 1) return false;

  const double kRRNoCap = 0.5 * std::numeric_limits<double>::max();
  const bool capRR = maxRR < kRRNoCap;
  for (int ihit = skipHitsFirst; ihit < n - skipHitsLast; ++ihit) {
    const AnaHitPD& h = part->Hits[2][ihit];
    if (capRR && h.ResidualRange > maxRR) continue;
    dedx.push_back(static_cast<double>(h.dEdx));
    rr.push_back(static_cast<double>(h.ResidualRange));
  }

  const bool dedxWindow = (dedxMinMeVcm > 0. && dedxMaxMeVcm > dedxMinMeVcm);
  if (dedxWindow) {
    std::vector<double> nd, nr;
    nd.reserve(dedx.size());
    nr.reserve(rr.size());
    for (size_t i = 0; i < dedx.size(); ++i) {
      if (dedx[i] >= dedxMinMeVcm && dedx[i] <= dedxMaxMeVcm) {
        nd.push_back(dedx[i]);
        nr.push_back(rr[i]);
      }
    }
    dedx.swap(nd);
    rr.swap(nr);
  }
  return static_cast<int>(dedx.size()) >= minPoints;
}

bool InteriorDedxRrSampleFromVectors(const std::vector<double>& rrIn, const std::vector<double>& dedxIn, int minPoints,
                                     std::vector<double>& dedx, std::vector<double>& rr, int skipFirst, int skipLast,
                                     double dedxMinMeVcm, double dedxMaxMeVcm) {
  dedx.clear();
  rr.clear();
  if (rrIn.size() != dedxIn.size() || rrIn.empty()) return false;
  if (skipFirst < 0) skipFirst = 0;
  if (skipLast < 0) skipLast = 0;

  const int n = static_cast<int>(rrIn.size());
  if (n < skipFirst + skipLast + 1) return false;

  for (int i = skipFirst; i < n - skipLast; ++i) {
    const double rrVal = rrIn[static_cast<size_t>(i)];
    const double dedxVal = dedxIn[static_cast<size_t>(i)];
    if (!std::isfinite(rrVal) || !std::isfinite(dedxVal) || !(rrVal > 0.)) continue;
    dedx.push_back(dedxVal);
    rr.push_back(rrVal);
  }

  const bool dedxWindow = (dedxMinMeVcm > 0. && dedxMaxMeVcm > dedxMinMeVcm);
  if (dedxWindow) {
    std::vector<double> nd, nr;
    nd.reserve(dedx.size());
    nr.reserve(rr.size());
    for (size_t i = 0; i < dedx.size(); ++i) {
      if (dedx[i] >= dedxMinMeVcm && dedx[i] <= dedxMaxMeVcm) {
        nd.push_back(dedx[i]);
        nr.push_back(rr[i]);
      }
    }
    dedx.swap(nd);
    rr.swap(nr);
  }
  return static_cast<int>(dedx.size()) >= minPoints;
}

double MeasuredTrackLengthCm(const AnaParticlePD* part, const std::vector<double>& rrInterior) {
  if (part && part->Length > 0.f && part->Length != -999.f) return static_cast<double>(part->Length);
  double mx = 0.;
  for (double r : rrInterior) {
    if (r > mx) mx = r;
  }
  return mx;
}

double ComputeScanL0FromInteriorRR(const std::vector<double>& rrInterior) {
  if (rrInterior.empty()) return 0.;
  double rrMin = std::numeric_limits<double>::infinity();
  for (double r : rrInterior) {
    if (!std::isfinite(r)) continue;
    rrMin = std::min(rrMin, r);
  }
  if (!std::isfinite(rrMin) || rrMin <= 0.) return 0.;
  return -rrMin;
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
                                                       bool computeMomentum, int skipHitsFirst, int skipHitsLast,
                                                       double dedxMinMeVcm, double dedxMaxMeVcm, double pdfPathCm,
                                                       bool computeDedxBiasDiagnostics = true,
                                                       double scanStepFineCm = 0., double lowPMomentumRefineGeV = 0.2) {
  pdAnaUtils::DEdxFreeRangeFitResult bad;
  if (!part || part->Hits[2].empty()) return bad;
  std::string ssparticle;
  Float_t mass = 0.f;
  if (!PdgToMassAndKeName(PDG, mass, ssparticle)) return bad;

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, maxRRForHits, minInteriorPoints, dedx, rr, skipHitsFirst, skipHitsLast,
                             dedxMinMeVcm, dedxMaxMeVcm))
    return bad;

  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return bad;

  TGraph* tg = new TGraph((int)dedx.size(), &rr[0], &dedx[0]);

  const double lenCm = MeasuredTrackLengthCm(part, rr);
  const double L0 = ComputeScanL0FromInteriorRR(rr);
  pdAnaUtils::DEdxFreeRangeFitResult result =
      pdAnaUtils::dEdxLikelihoodFreeRangeFit(tg, tg_ke, mass, L0, Lmax, step, lenCm, computeMomentum, pdfPathCm,
                                             computeDedxBiasDiagnostics, scanStepFineCm, lowPMomentumRefineGeV);

  delete tg;
  return result;
}

bool BuildPionFreeRangeLogLikelihoodVsMomentumCurveInternal(AnaParticlePD* part, double Lmax, double step,
                                                            int minInteriorPoints, int skipHitsFirst, int skipHitsLast,
                                                            double dedxMinMeVcm, double dedxMaxMeVcm, double pdfPathCm,
                                                            std::vector<double>& pGeV, std::vector<double>& logL,
                                                            double scanStepFineCm, double lowPMomentumRefineGeV) {
  pGeV.clear();
  logL.clear();
  if (!part || part->Hits[2].empty()) return false;
  constexpr Int_t kPdg = 211;
  std::string ssparticle;
  Float_t mass = 0.f;
  if (!PdgToMassAndKeName(kPdg, mass, ssparticle)) return false;

  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), minInteriorPoints, dedx, rr, skipHitsFirst,
                           skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm))
    return false;

  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return false;

  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());
  const double lenCm = MeasuredTrackLengthCm(part, rr);

  if (!tg || tg->GetN() < 1 || step <= 0.) {
    delete tg;
    return false;
  }
  if (!std::isfinite(pdfPathCm) || pdfPathCm <= 0.) {
    delete tg;
    return false;
  }

  const double L0 = ComputeScanL0FromInteriorRR(rr);
  if (!(Lmax >= L0)) {
    delete tg;
    return false;
  }

  const double width = pdfPathCm;
  TF1* pdf = FreeRangeLikelihoodPdf();

  const auto scanAtStep = [&](double stepUse, std::vector<std::pair<double, double>>& pairOut) -> bool {
    pairOut.clear();
    pairOut.reserve(static_cast<size_t>(std::max(1.0, (Lmax - L0) / stepUse) + 2.));
    for (double L = L0; L <= Lmax + 1e-9; L += stepUse) {
      double ll = 0.;
      for (int i = 0; i < tg->GetN(); i++) {
        const double range = tg->GetPointX(i) + L;
        const double dEdx = tg->GetPointY(i);
        if (!(range > 0.) || !std::isfinite(range)) continue;
        const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
        if (ke < 0. || !std::isfinite(ke)) continue;
        double par[5] = {0., 0., 0., 0., 0.};
        if (!BuildDedxPdfParams(ke, mass, width, par)) continue;
        pdf->SetParameters(par);
        const double pval = pdf->Eval(dEdx);
        if (pval == 0.) continue;
        ll += std::log(pval);
      }
      if (!std::isfinite(ll)) continue;
      const double R_eff = lenCm + L;
      if (!(R_eff > 0.) || !std::isfinite(R_eff)) continue;
      const double p = pdMomShared::RangeCmToMomentumGeV(R_eff, kPdg, tg_ke, mass);
      if (!std::isfinite(p) || p <= 0.) continue;
      pairOut.push_back({p, ll});
    }
    return pairOut.size() >= 2u;
  };

  std::vector<std::pair<double, double>> pairs;
  if (!scanAtStep(step, pairs)) {
    delete tg;
    return false;
  }
  bool needFine = false;
  if (lowPMomentumRefineGeV > 0.) {
    for (const auto& pr : pairs) {
      if (pr.first < lowPMomentumRefineGeV) {
        needFine = true;
        break;
      }
    }
  }
  if (needFine && scanStepFineCm > 0. && scanStepFineCm + 1e-12 < step) {
    std::vector<std::pair<double, double>> fine;
    if (scanAtStep(scanStepFineCm, fine)) pairs.swap(fine);
  }
  delete tg;

  if (pairs.size() < 2u) return false;

  std::sort(pairs.begin(), pairs.end(),
            [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first < b.first; });

  for (const auto& pr : pairs) {
    if (!pGeV.empty() && std::abs(pr.first - pGeV.back()) <= 1e-12 * std::max(1.0, std::abs(pr.first))) {
      if (pr.second > logL.back()) logL.back() = pr.second;
    } else {
      pGeV.push_back(pr.first);
      logL.push_back(pr.second);
    }
  }
  return pGeV.size() >= 2u;
}

} // namespace

//***************************************************************
bool pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(AnaParticlePD* part, double Lmax, double step,
                                                                int minInteriorPoints, int skipHitsFirst,
                                                                int skipHitsLast, double dedxMinMeVcm,
                                                                double dedxMaxMeVcm, double pdfPathCm,
                                                                std::vector<double>& pGeV, std::vector<double>& logL,
                                                                double scanStepFineCm, double lowPMomentumRefineGeV) {
//***************************************************************
  return BuildPionFreeRangeLogLikelihoodVsMomentumCurveInternal(part, Lmax, step, minInteriorPoints, skipHitsFirst,
                                                                skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm, pdfPathCm,
                                                                pGeV, logL, scanStepFineCm, lowPMomentumRefineGeV);
}

//***************************************************************
bool pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurveFromVectors(
    const std::vector<double>& rrCm, const std::vector<double>& dedxMeVcm, double measuredTrackLengthCm, double Lmax,
    double step, int minInteriorPoints, int skipFirst, int skipLast, double dedxMinMeVcm, double dedxMaxMeVcm,
    double pdfPathCm, std::vector<double>& pGeV, std::vector<double>& logL, double scanStepFineCm,
    double lowPMomentumRefineGeV) {
//***************************************************************
  pGeV.clear();
  logL.clear();
  constexpr Int_t kPdg = 211;
  std::string ssparticle;
  Float_t mass = 0.f;
  if (!PdgToMassAndKeName(kPdg, mass, ssparticle)) return false;

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSampleFromVectors(rrCm, dedxMeVcm, minInteriorPoints, dedx, rr, skipFirst, skipLast, dedxMinMeVcm,
                                       dedxMaxMeVcm))
    return false;

  TGraph* tg_ke = KeVsRangeGraphCached(ssparticle);
  if (!tg_ke) return false;

  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());
  const double lenCm = (std::isfinite(measuredTrackLengthCm) && measuredTrackLengthCm > 0.)
                           ? measuredTrackLengthCm
                           : MeasuredTrackLengthCm(nullptr, rr);

  if (!tg || tg->GetN() < 1 || step <= 0.) {
    delete tg;
    return false;
  }
  if (!std::isfinite(pdfPathCm) || pdfPathCm <= 0.) {
    delete tg;
    return false;
  }

  const double L0 = ComputeScanL0FromInteriorRR(rr);
  if (!(Lmax >= L0)) {
    delete tg;
    return false;
  }

  const double width = pdfPathCm;
  TF1* pdf = FreeRangeLikelihoodPdf();

  const auto scanAtStep = [&](double stepUse, std::vector<std::pair<double, double>>& pairOut) -> bool {
    pairOut.clear();
    pairOut.reserve(static_cast<size_t>(std::max(1.0, (Lmax - L0) / stepUse) + 2.));
    for (double L = L0; L <= Lmax + 1e-9; L += stepUse) {
      double ll = 0.;
      for (int i = 0; i < tg->GetN(); i++) {
        const double range = tg->GetPointX(i) + L;
        const double dEdx = tg->GetPointY(i);
        if (!(range > 0.) || !std::isfinite(range)) continue;
        const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
        if (ke < 0. || !std::isfinite(ke)) continue;
        double par[5] = {0., 0., 0., 0., 0.};
        if (!BuildDedxPdfParams(ke, mass, width, par)) continue;
        pdf->SetParameters(par);
        const double pval = pdf->Eval(dEdx);
        if (pval == 0.) continue;
        ll += std::log(pval);
      }
      if (!std::isfinite(ll)) continue;
      const double R_eff = lenCm + L;
      if (!(R_eff > 0.) || !std::isfinite(R_eff)) continue;
      const double p = pdMomShared::RangeCmToMomentumGeV(R_eff, kPdg, tg_ke, mass);
      if (!std::isfinite(p) || p <= 0.) continue;
      pairOut.push_back({p, ll});
    }
    return pairOut.size() >= 2u;
  };

  std::vector<std::pair<double, double>> pairs;
  if (!scanAtStep(step, pairs)) {
    delete tg;
    return false;
  }
  bool needFine = false;
  if (lowPMomentumRefineGeV > 0.) {
    for (const auto& pr : pairs) {
      if (pr.first < lowPMomentumRefineGeV) {
        needFine = true;
        break;
      }
    }
  }
  if (needFine && scanStepFineCm > 0. && scanStepFineCm + 1e-12 < step) {
    std::vector<std::pair<double, double>> fine;
    if (scanAtStep(scanStepFineCm, fine)) pairs.swap(fine);
  }
  delete tg;

  if (pairs.size() < 2u) return false;

  std::sort(pairs.begin(), pairs.end(),
            [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first < b.first; });

  for (const auto& pr : pairs) {
    if (!pGeV.empty() && std::abs(pr.first - pGeV.back()) <= 1e-12 * std::max(1.0, std::abs(pr.first))) {
      if (pr.second > logL.back()) logL.back() = pr.second;
    } else {
      pGeV.push_back(pr.first);
      logL.push_back(pr.second);
    }
  }
  return pGeV.size() >= 2u;
}

static TF1* FreeRangeLikelihoodPdf() {
  static TF1* pdf =
      new TF1("pdf_dedx_freerange_global", pdAnaUtils::dEdxPDF, -10., 20., 5);
  return pdf;
}

static bool BuildDedxPdfParams(double ke, Float_t mass, double width, double* par) {
  if (!par || !std::isfinite(ke) || ke <= 0. || !std::isfinite(static_cast<double>(mass)) || mass <= 0.) return false;
  const double gamma = (ke / mass) + 1.0;
  if (!std::isfinite(gamma) || gamma <= 1.0) return false;
  const double beta2 = 1.0 - (1.0 / (gamma * gamma));
  if (!std::isfinite(beta2) || beta2 <= 0.) return false;
  const double beta = std::sqrt(beta2);
  const double xi = pdAnaUtils::GetLandauXi(ke, width, mass);
  const double Wmax = pdAnaUtils::GetWmax(ke, mass);
  if (!std::isfinite(xi) || !std::isfinite(Wmax) || xi <= 0. || Wmax <= 0.) return false;
  const double kappa = xi / Wmax;
  const double dEdx_BB = pdAnaUtils::GetdEdxBetheBloch(ke, mass);
  if (!std::isfinite(kappa) || !std::isfinite(dEdx_BB) || dEdx_BB <= 0.) return false;
  par[0] = kappa;
  par[1] = beta * beta;
  par[2] = xi;
  par[3] = dEdx_BB;
  par[4] = width;
  return true;
}

static double ExpectedDedxFromLikelihoodPdfMode(double ke, Float_t mass, double width) {
  double par[5] = {0., 0., 0., 0., 0.};
  if (!BuildDedxPdfParams(ke, mass, width, par)) return -1.;

  const double dEdxBB = par[3];
  const double xmin = std::max(0.02, dEdxBB - 8.0);
  const double xmax = std::min(60.0, dEdxBB + 12.0);
  if (!(xmax > xmin)) return -1.;

  auto pdfEval = [&](double x) {
    double xx[1] = {x};
    return pdAnaUtils::dEdxPDF(xx, par);
  };

  const int nScan = 400;
  const double step = (xmax - xmin) / static_cast<double>(nScan);
  double xBest = xmin;
  double pBest = -1.;
  for (int i = 0; i <= nScan; ++i) {
    const double x = xmin + step * static_cast<double>(i);
    const double p = pdfEval(x);
    if (std::isfinite(p) && p > pBest) {
      pBest = p;
      xBest = x;
    }
  }
  if (!std::isfinite(pBest) || pBest <= 0.) return -1.;

  double halfWindow = std::max(0.02, 3.0 * step);
  for (int pass = 0; pass < 2; ++pass) {
    const double lo = std::max(xmin, xBest - halfWindow);
    const double hi = std::min(xmax, xBest + halfWindow);
    const int nFine = 200;
    const double fineStep = (hi - lo) / static_cast<double>(nFine);
    for (int i = 0; i <= nFine; ++i) {
      const double x = lo + fineStep * static_cast<double>(i);
      const double p = pdfEval(x);
      if (std::isfinite(p) && p > pBest) {
        pBest = p;
        xBest = x;
      }
    }
    halfWindow *= 0.35;
  }

  return xBest;
}

static pdAnaUtils::DEdxFreeRangeFitResult RunFreeRangeScan(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0, double Lmax,
                                                         double step, double measuredTrackLengthCm, bool computeMomentum,
                                                         double pdfPathCm, bool computeDedxBiasDiagnostics,
                                                         double scanStepFineCm, double lowPMomentumRefineGeV) {
  pdAnaUtils::DEdxFreeRangeFitResult out;
  if (!tg || !tg_ke || tg->GetN() < 1 || step <= 0 || Lmax < L0) return out;
  if (!std::isfinite(pdfPathCm) || pdfPathCm <= 0.) return out;

  const double width = pdfPathCm;
  TF1* pdf = FreeRangeLikelihoodPdf();
  std::vector<double> L_v;
  std::vector<double> Likelihood_v;

  const auto fillGrid = [&](double stepUse, std::vector<double>& Lout, std::vector<double>& LLout, bool* sawLowP) {
    Lout.clear();
    LLout.clear();
    if (sawLowP) *sawLowP = false;
    for (double L = L0; L <= Lmax + 1e-9; L += stepUse) {
      double likelihood = 0.;
      for (int i = 0; i < tg->GetN(); i++) {
        const double range = tg->GetPointX(i) + L;
        const double dEdx = tg->GetPointY(i);
        if (!(range > 0.) || !std::isfinite(range)) continue;
        const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
        if (ke < 0. || !std::isfinite(ke)) continue;
        double par[5] = {0., 0., 0., 0., 0.};
        if (!BuildDedxPdfParams(ke, mass, width, par)) continue;
        pdf->SetParameters(par);
        const double pval = pdf->Eval(dEdx);
        if (pval == 0.) continue;
        likelihood += std::log(pval);
      }
      Lout.push_back(L);
      LLout.push_back(likelihood);
      if (sawLowP && computeMomentum && measuredTrackLengthCm > 0. && std::isfinite(measuredTrackLengthCm) &&
          lowPMomentumRefineGeV > 0.) {
        const double R_eff = measuredTrackLengthCm + L;
        const double pTry = pdMomShared::RangeCmToMomentumGeV(R_eff, 211, tg_ke, mass);
        if (std::isfinite(pTry) && pTry > 0. && pTry < lowPMomentumRefineGeV) *sawLowP = true;
      }
    }
  };

  bool probeLowP = false;
  fillGrid(step, L_v, Likelihood_v, &probeLowP);
  if (scanStepFineCm > 0. && scanStepFineCm + 1e-12 < step && probeLowP && lowPMomentumRefineGeV > 0.) {
    std::vector<double> Lfine, LLfine;
    bool dummy = false;
    fillGrid(scanStepFineCm, Lfine, LLfine, &dummy);
    if (!LLfine.empty()) {
      L_v.swap(Lfine);
      Likelihood_v.swap(LLfine);
    }
  }

  if (Likelihood_v.empty()) return out;

  auto it = std::max_element(Likelihood_v.begin(), Likelihood_v.end());
  const int index = static_cast<int>(std::distance(Likelihood_v.begin(), it));
  out.logLikelihood = static_cast<Float_t>(Likelihood_v[index]);
  out.bestOffsetCm = static_cast<Float_t>(L_v[index]);

  if (computeDedxBiasDiagnostics) {
    // Compute dEdx bias distribution and fit Gaussian at best offset
    std::vector<double> dedx_bias_values;
    const double L_best = L_v[index];
    for (int i = 0; i < tg->GetN(); i++) {
      const double range = tg->GetPointX(i) + L_best;
      const double dEdx = tg->GetPointY(i);
      if (!(range > 0.) || !std::isfinite(range)) continue;
      const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, range);
      if (ke < 0. || !std::isfinite(ke)) continue;
      const double dEdx_expected = ExpectedDedxFromLikelihoodPdfMode(ke, mass, width);
      if (!std::isfinite(dEdx_expected) || dEdx_expected <= 0.) continue;
      dedx_bias_values.push_back(dEdx - dEdx_expected);
    }

    // Fit Gaussian to dEdx bias distribution in a robust central window
    if (dedx_bias_values.size() >= 3) {
      // 0.1 bin width in [-1, 1]
      TH1F hist_bias("h_dedx_bias", "dEdx bias distribution", 20, -1., 1.);
      for (double bias : dedx_bias_values) {
        hist_bias.Fill(bias);
      }

      if (hist_bias.GetEntries() >= 3 && std::isfinite(hist_bias.GetRMS()) && hist_bias.GetRMS() > 0.) {
        TF1* gaus = new TF1("gaus_fit", "gaus", -1., 1.);
        gaus->SetParameter(0, std::max(1.0, static_cast<double>(hist_bias.GetMaximum())));
        gaus->SetParameter(1, hist_bias.GetMean());
        gaus->SetParameter(2, std::max(0.05, static_cast<double>(hist_bias.GetRMS())));
        gaus->SetParLimits(1, -1., 1.);
        gaus->SetParLimits(2, 0.02, 1.0);

        const int fitStatus = hist_bias.Fit(gaus, "RQ0", "", -1., 1.);
        if (fitStatus == 0 && std::isfinite(gaus->GetParameter(1)) && std::isfinite(gaus->GetParameter(2)) &&
            gaus->GetParameter(2) > 0.) {
          out.meanDedxBias = static_cast<Float_t>(gaus->GetParameter(1));
          out.sigmaDedxBias = static_cast<Float_t>(gaus->GetParameter(2));
          out.dedxFitOk = 1;
        } else {
          out.meanDedxBias = static_cast<Float_t>(hist_bias.GetMean());
          out.sigmaDedxBias = static_cast<Float_t>(hist_bias.GetRMS());
          out.dedxFitOk = 0;
        }
        delete gaus;
      } else {
        out.meanDedxBias = static_cast<Float_t>(hist_bias.GetMean());
        out.sigmaDedxBias = static_cast<Float_t>(hist_bias.GetRMS());
        out.dedxFitOk = 0;
      }
    }
  }

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
Float_t pdAnaUtils::GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG, int skipHitsFirst, int skipHitsLast,
                                      double dedxMinMeVcm, double dedxMaxMeVcm) {
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
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), 2, dedx, rr, skipHitsFirst, skipHitsLast,
                             dedxMinMeVcm, dedxMaxMeVcm)) {
    return -999.;
  }
  TGraph* tg = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;

  return result;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR, int skipHitsFirst,
                                             int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm) {
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
  if (!InteriorDedxRrSample(part, maxRR, 2, dedx, rr, skipHitsFirst, skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm)) {
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
                                                            bool computeMomentum, double pdfPathCm,
                                                            bool computeDedxBiasDiagnostics, double scanStepFineCm,
                                                            double lowPMomentumRefineGeV) {
//***************************************************************
  return RunFreeRangeScan(tg, tg_ke, mass, L0, Lmax, step, measuredTrackLengthCm, computeMomentum, pdfPathCm,
                          computeDedxBiasDiagnostics, scanStepFineCm, lowPMomentumRefineGeV);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke,
                                                                Float_t mass){
//***************************************************************
  const pdAnaUtils::DEdxFreeRangeFitResult r =
      dEdxLikelihoodFreeRangeFit(tg, tg_ke, mass, 0., 10., 0.1, -1., false, 0.65);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
std::pair<Float_t, Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG, int skipHitsFirst,
                                                                   int skipHitsLast, double dedxMinMeVcm,
                                                                   double dedxMaxMeVcm, double pdfPathCm) {
//***************************************************************

  const pdAnaUtils::DEdxFreeRangeFitResult r = ParticleFreeRangeFit(
      part, PDG, 10., 0.1, std::numeric_limits<double>::max(), 1, false, skipHitsFirst, skipHitsLast, dedxMinMeVcm,
      dedxMaxMeVcm, pdfPathCm, false);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
pdAnaUtils::DEdxFreeRangeFitResult pdAnaUtils::GetdEdxLikelihoodFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax,
                                                                             double step, int minInteriorPoints,
                                                                             int skipHitsFirst, int skipHitsLast,
                                                                             double dedxMinMeVcm, double dedxMaxMeVcm,
                                                                             double pdfPathCm,
                                                                             bool computeDedxBiasDiagnostics,
                                                                             double scanStepFineCm,
                                                                             double lowPMomentumRefineGeV) {
//***************************************************************
  return ParticleFreeRangeFit(part, PDG, Lmax, step, std::numeric_limits<double>::max(), minInteriorPoints, true,
                              skipHitsFirst, skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm, pdfPathCm,
                              computeDedxBiasDiagnostics, scanStepFineCm, lowPMomentumRefineGeV);
}

//***************************************************************
std::pair<Float_t, Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG,
                                                                          const double maxRR, int skipHitsFirst,
                                                                          int skipHitsLast, double dedxMinMeVcm,
                                                                          double dedxMaxMeVcm, double pdfPathCm) {
//***************************************************************

  const pdAnaUtils::DEdxFreeRangeFitResult r =
      ParticleFreeRangeFit(part, PDG, 10., 0.1, maxRR, 1, false, skipHitsFirst, skipHitsLast, dedxMinMeVcm,
                           dedxMaxMeVcm, pdfPathCm, false);
  return std::make_pair(r.logLikelihood, r.bestOffsetCm);
}

//***************************************************************
pdAnaUtils::DEdxFreeRangeFitResult pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG,
                                                                                     const double maxRR, double Lmax,
                                                                                     double step, int minInteriorPoints,
                                                                                     int skipHitsFirst, int skipHitsLast,
                                                                                     double dedxMinMeVcm,
                                                                                     double dedxMaxMeVcm,
                                                                                     double pdfPathCm,
                                                                                     bool computeDedxBiasDiagnostics,
                                                                                     double scanStepFineCm,
                                                                                     double lowPMomentumRefineGeV) {
//***************************************************************
  return ParticleFreeRangeFit(part, PDG, Lmax, step, maxRR, minInteriorPoints, true, skipHitsFirst, skipHitsLast,
                              dedxMinMeVcm, dedxMaxMeVcm, pdfPathCm, computeDedxBiasDiagnostics, scanStepFineCm,
                              lowPMomentumRefineGeV);
}

//***************************************************************
TMultiGraph* pdAnaUtils::MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax, double step,
                                                             int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                                             double dedxMaxMeVcm, int minInteriorPoints,
                                                             double pdfPathCm, const char* xAxisTitle) {
//***************************************************************
  if (!part || part->Hits[2].empty()) return nullptr;
  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), minInteriorPoints, dedx, rr, skipHitsFirst,
                            skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm))
    return nullptr;

  TGraph* tg_ke = KeVsRangeGraphCached("pion");
  if (!tg_ke) return nullptr;
  TGraph* tgFit = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());
  const double lenCm = MeasuredTrackLengthCm(part, rr);
  const double L0 = ComputeScanL0FromInteriorRR(rr);
  const DEdxFreeRangeFitResult fit = pdAnaUtils::dEdxLikelihoodFreeRangeFit(
      tgFit, tg_ke, 139.57f, L0, Lmax, step, lenCm, true, pdfPathCm);
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

  constexpr double kRefCurveStepCm = 1.0;
  constexpr double kLowRStartCm = 0.02;
  constexpr double kPionMassMeV = 139.57;
  const double kPdfWidth = pdfPathCm;
  TGraph* tg_ke_pion = KeVsRangeGraphCached("pion");
  if (!tg_ke_pion) {
    delete gMeas;
    delete gShift;
    return nullptr;
  }

  double rrMaxMeasured = 0.;
  for (double r : rr) {
    if (std::isfinite(r)) rrMaxMeasured = std::max(rrMaxMeasured, r);
  }
  double rrMaxShifted = 0.;
  for (double r : rrShift) {
    if (std::isfinite(r)) rrMaxShifted = std::max(rrMaxShifted, r);
  }
  const double rrMaxFromScan = rrMaxMeasured + std::max(0.0, Lmax);
  const double kRefCurveRrMaxCm = std::max({kLowRStartCm, rrMaxMeasured, rrMaxShifted, rrMaxFromScan});

  std::vector<double> tplX, tplY;
  for (double rrRef = kLowRStartCm; rrRef <= kRefCurveRrMaxCm + 1e-9; rrRef += kRefCurveStepCm) {
    const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke_pion, rrRef);
    if (ke < 0. || !std::isfinite(ke)) continue;
    const double yexp = ExpectedDedxFromLikelihoodPdfMode(ke, static_cast<Float_t>(kPionMassMeV), kPdfWidth);
    if (!std::isfinite(yexp) || yexp <= 0.) continue;
    tplX.push_back(rrRef);
    tplY.push_back(yexp);
  }
  if (tplX.empty()) {
    delete gMeas;
    delete gShift;
    return nullptr;
  }

  auto* gTpl = new TGraph(static_cast<int>(tplX.size()), tplX.data(), tplY.data());
  gTpl->SetName("g_pion_dedx_pdf_expected");
  gTpl->SetTitle("pion expected dE/dx from free-range likelihood PDF mode");
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
TH1F* pdAnaUtils::MakePionFreeRangeDedxBiasHistogram(AnaParticlePD* part, double Lmax, double step, int skipHitsFirst,
                                                     int skipHitsLast, double dedxMinMeVcm, double dedxMaxMeVcm,
                                                     int minInteriorPoints, double pdfPathCm, const char* histTitle) {
//***************************************************************
  if (!part || part->Hits[2].empty()) return nullptr;
  if (CollectionPlaneResidualRangeLooksUnset(part)) pdAnaUtils::ComputeResidualRange(part);

  std::vector<double> dedx, rr;
  if (!InteriorDedxRrSample(part, std::numeric_limits<double>::max(), minInteriorPoints, dedx, rr, skipHitsFirst,
                             skipHitsLast, dedxMinMeVcm, dedxMaxMeVcm))
    return nullptr;

  TGraph* tg_ke = KeVsRangeGraphCached("pion");
  if (!tg_ke) return nullptr;

  TGraph* tgFit = new TGraph(static_cast<int>(dedx.size()), rr.data(), dedx.data());
  const double lenCm = MeasuredTrackLengthCm(part, rr);
  const double L0 = ComputeScanL0FromInteriorRR(rr);
  const DEdxFreeRangeFitResult fit = pdAnaUtils::dEdxLikelihoodFreeRangeFit(
      tgFit, tg_ke, 139.57f, L0, Lmax, step, lenCm, true, pdfPathCm);
  delete tgFit;
  if (!std::isfinite(static_cast<double>(fit.bestOffsetCm)) || fit.bestOffsetCm <= -998.f) return nullptr;

  const double L = static_cast<double>(fit.bestOffsetCm);
  std::vector<double> bias;
  bias.reserve(dedx.size());
  for (size_t i = 0; i < dedx.size(); ++i) {
    const double rrShifted = rr[i] + L;
    if (!(rrShifted > 0.) || !std::isfinite(rrShifted)) continue;
    const double ke = pdAnaUtils::KineticEnergyMeVFromResidualRangeCm(tg_ke, rrShifted);
    if (ke < 0. || !std::isfinite(ke)) continue;
    const double dedxExpected = ExpectedDedxFromLikelihoodPdfMode(ke, 139.57f, static_cast<Float_t>(pdfPathCm));
    if (!std::isfinite(dedxExpected) || dedxExpected <= 0.) continue;
    bias.push_back(dedx[i] - dedxExpected);
  }
  if (bias.size() < 3u) return nullptr;

  auto* h = new TH1F("h_pion_freerange_dedx_bias",
                     (histTitle && histTitle[0]) ? histTitle : "#Delta(dE/dx)=measured-expected(PDF mode);#Delta(dE/dx) [MeV/cm];Entries",
                     20, -1., 1.);
  h->SetDirectory(nullptr);
  for (double v : bias) h->Fill(v);

  if (h->GetEntries() >= 3 && std::isfinite(h->GetRMS()) && h->GetRMS() > 0.) {
    TF1* gaus = new TF1("gaus_dedx_bias_hist", "gaus", -1., 1.);
    gaus->SetParameter(0, std::max(1.0, static_cast<double>(h->GetMaximum())));
    gaus->SetParameter(1, h->GetMean());
    gaus->SetParameter(2, std::max(0.05, static_cast<double>(h->GetRMS())));
    gaus->SetParLimits(1, -1., 1.);
    gaus->SetParLimits(2, 0.02, 1.0);
    h->Fit(gaus, "RQ0", "", -1., 1.);
    delete gaus;
  }

  return h;
}
