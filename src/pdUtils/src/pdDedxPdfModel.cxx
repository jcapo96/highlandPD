#include "pdDedxPdfModel.hxx"

#include <Math/VavilovAccurate.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>

namespace pdDedxPdf {

namespace {

Model gModel;

// Physical constants, as in lambda_feasibility.py (PDG; Sternheimer density-effect parameters for liquid argon).
constexpr double kK = 0.307075;  // MeV cm2/mol
constexpr double kZoverA = 18.0 / 39.948;
constexpr double kI = 188.0e-6;  // MeV
constexpr double kMe = 0.51099895;
constexpr double kSX0 = 0.2, kSX1 = 3.0, kSA = 0.19559, kSK = 3.0, kSC = 5.2146;

// PDF tables: cells on a grid of ln(KE) and ln(pitch), each holding the PDF on the dE/dx window.
constexpr double kU0 = -0.6931471805599453;  // ln(0.5 MeV)
constexpr double kDU = 0.01;
constexpr int kNU = 1064;                    // up to ~20 GeV
constexpr double kV0 = -1.3862943611198906;  // ln(0.25 cm)
constexpr double kDV = 0.02;
constexpr int kNV = 185;                     // up to ~10 cm
constexpr double kDX = 0.005;                // MeV/cm, dE/dx step of the tables
constexpr double kMinPdf = 1e-12;

struct Cache {
  std::map<double, std::unordered_map<long, std::vector<float>>> cells;  // by mass
  std::map<double, std::pair<std::vector<double>, std::vector<double>>> range;  // by mass: (KE, R)
};
Cache& cache() {
  static Cache c;
  return c;
}

double DensityCorrection(double bg) {
  const double x = std::log10(bg);
  if (x >= kSX1) return 2. * std::log(10.) * x - kSC;
  if (x >= kSX0) return 2. * std::log(10.) * x - kSC + kSA * std::pow(kSX1 - x, kSK);
  return 0.;
}

struct Kin {
  double beta2, gamma, wmax;
};
Kin Kinematics(double ke, double mass) {
  Kin k;
  k.gamma = 1. + ke / mass;
  k.beta2 = 1. - 1. / (k.gamma * k.gamma);
  const double bg2 = k.beta2 * k.gamma * k.gamma;
  k.wmax = 2. * kMe * bg2 / (1. + 2. * k.gamma * kMe / mass + (kMe / mass) * (kMe / mass));
  return k;
}

const std::pair<std::vector<double>, std::vector<double>>& RangeTable(double mass) {
  auto& tab = cache().range[mass];
  if (!tab.first.empty()) return tab;
  constexpr int n = 4000;
  const double lo = std::log10(0.1), hi = std::log10(5000.);
  std::vector<double> t(n), r(n, 0.);
  for (int i = 0; i < n; ++i) t[i] = std::pow(10., lo + (hi - lo) * i / (n - 1));
  for (int i = 1; i < n; ++i) r[i] = r[i - 1] + (t[i] - t[i - 1]) / MeanDedx(0.5 * (t[i] + t[i - 1]), mass);
  tab.first = std::move(t);
  tab.second = std::move(r);
  return tab;
}

// Energy-loss PDF (no detector resolution) at dE/dx x, in the Landau/Vavilov/Gauss regimes of the legacy dEdxPDF,
// with the Gaussian branch normalised.
struct Straggling {
  double kappa, beta2, a, b;  // dE/dx = b + a * lambda
  bool landau, gauss;
  double gmu, gsigma;
  ROOT::Math::VavilovAccurate* vav = nullptr;
  double operator()(double x) const {
    const double y = (x - b) / a;
    if (landau) return TMath::Landau(y) / a;
    if (gauss) return TMath::Gaus(y, gmu, gsigma, true) / a;
    return vav->Pdf(y) / a;
  }
};

std::vector<float> BuildCell(int iu, int iv, double mass) {
  const Model& m = gModel;
  const double ke = std::exp(kU0 + kDU * iu);
  const double pitch = std::exp(kV0 + kDV * iv);
  const Kin k = Kinematics(ke, mass);
  const double mean = MeanDedx(ke, mass);
  const double xi = 0.5 * kK * kZoverA * m.densityGcm3 * pitch / k.beta2;
  static ROOT::Math::VavilovAccurate vav;
  Straggling s;
  s.kappa = xi / k.wmax;
  s.beta2 = k.beta2;
  s.a = xi / pitch;
  s.b = mean + s.a * (0.422784 + k.beta2 + std::log(s.kappa));
  s.landau = s.kappa < 0.01;
  s.gauss = s.kappa > 10.;
  if (s.gauss) {
    s.gmu = ROOT::Math::Vavilov::Mean(s.kappa, s.beta2);
    s.gsigma = std::sqrt(ROOT::Math::Vavilov::Variance(s.kappa, s.beta2));
  } else if (!s.landau) {
    vav.SetKappaBeta2(s.kappa, s.beta2);
    s.vav = &vav;
  }

  const double sigma = std::sqrt(std::pow(m.resolutionFrac * mean, 2) + m.resolutionAbs * m.resolutionAbs);
  const int nw = static_cast<int>(std::lround((m.dedxMax - m.dedxMin) / kDX)) + 1;
  std::vector<double> g(nw, 0.);
  if (sigma < 0.5 * kDX) {
    for (int j = 0; j < nw; ++j) g[j] = s(m.dedxMin + kDX * j);
  } else {
    // Convolution on a grid that extends 6 sigma beyond the window on both sides.
    const int ext = static_cast<int>(std::ceil(6. * sigma / kDX));
    const int nf = nw + 2 * ext;
    std::vector<double> f(nf);
    for (int j = 0; j < nf; ++j) {
      const double x = m.dedxMin + kDX * (j - ext);
      f[j] = (x > 0.) ? s(x) : 0.;
    }
    std::vector<double> kern(2 * ext + 1);
    for (int d = -ext; d <= ext; ++d) kern[d + ext] = TMath::Gaus(kDX * d, 0., sigma, true) * kDX;
    for (int j = 0; j < nw; ++j) {
      double acc = 0.;
      for (int d = -ext; d <= ext; ++d) acc += f[j + ext - d] * kern[d + ext];
      g[j] = acc;
    }
  }
  if (m.truncate) {
    double norm = 0.;
    for (int j = 0; j < nw; ++j) norm += g[j] * ((j == 0 || j == nw - 1) ? 0.5 : 1.) * kDX;
    if (norm > 1e-300)
      for (double& v : g) v /= norm;
  }
  const double unif = (m.dedxMax > m.dedxMin) ? 1. / (m.dedxMax - m.dedxMin) : 0.;
  std::vector<float> out(nw);
  for (int j = 0; j < nw; ++j) out[j] = static_cast<float>((1. - m.outlierFrac) * g[j] + m.outlierFrac * unif);
  return out;
}

const std::vector<float>& Cell(int iu, int iv, double mass) {
  auto& cells = cache().cells[mass];
  const long key = static_cast<long>(iu) * 1000L + iv;
  auto it = cells.find(key);
  if (it != cells.end()) return it->second;
  return cells.emplace(key, BuildCell(iu, iv, mass)).first->second;
}

double CellValue(const std::vector<float>& c, double x) {
  const double fx = (x - gModel.dedxMin) / kDX;
  if (fx <= 0.) return c.front();
  const int j = static_cast<int>(fx);
  if (j >= static_cast<int>(c.size()) - 1) return c.back();
  const double w = fx - j;
  return (1. - w) * c[j] + w * c[j + 1];
}

}  // namespace

void SetModel(const Model& m) {
  gModel = m;
  cache().cells.clear();
  cache().range.clear();
}
const Model& GetModel() { return gModel; }
bool Enabled() { return gModel.enabled; }

double MeanDedx(double keMeV, double massMeV) {
  const Kin k = Kinematics(keMeV, massMeV);
  const double bg = std::sqrt(k.beta2) * k.gamma;
  const double arg = 2. * kMe * bg * bg * k.wmax / (kI * kI);
  return kK * kZoverA * gModel.densityGcm3 / k.beta2 * (0.5 * std::log(arg) - k.beta2 - DensityCorrection(bg) / 2.);
}

double CsdaRangeCm(double keMeV, double massMeV) {
  const auto& tab = RangeTable(massMeV);
  const auto& t = tab.first;
  const auto& r = tab.second;
  if (!(keMeV > t.front())) return 0.;
  const size_t n = t.size();
  if (keMeV >= t.back()) return r[n - 1] + (keMeV - t[n - 1]) * (r[n - 1] - r[n - 2]) / (t[n - 1] - t[n - 2]);
  const size_t i = std::upper_bound(t.begin(), t.end(), keMeV) - t.begin();
  return r[i - 1] + (keMeV - t[i - 1]) * (r[i] - r[i - 1]) / (t[i] - t[i - 1]);
}

double KeFromRangeCm(double rangeCm, double massMeV) {
  const auto& tab = RangeTable(massMeV);
  const auto& t = tab.first;
  const auto& r = tab.second;
  if (!(rangeCm > 0.)) return t.front();
  const size_t n = r.size();
  if (rangeCm >= r.back()) return t[n - 1] + (rangeCm - r[n - 1]) * (t[n - 1] - t[n - 2]) / (r[n - 1] - r[n - 2]);
  const size_t i = std::upper_bound(r.begin(), r.end(), rangeCm) - r.begin();
  return t[i - 1] + (rangeCm - r[i - 1]) * (t[i] - t[i - 1]) / (r[i] - r[i - 1]);
}

double MomentumGeVFromRangeCm(double rangeCm, double massMeV) {
  if (!(rangeCm > 0.) || !std::isfinite(rangeCm)) return -999.;
  const double ke = KeFromRangeCm(rangeCm, massMeV);
  return std::sqrt(ke * ke + 2. * massMeV * ke) / 1000.;
}

double LogPdf(double dedx, double keMeV, double massMeV, double pitchCm) {
  const Model& m = gModel;
  if (!m.perHitPitch || !(pitchCm > 0.) || !std::isfinite(pitchCm)) pitchCm = m.fixedPitchCm;
  const double fu = std::min(std::max((std::log(std::max(keMeV, 1e-3)) - kU0) / kDU, 0.), kNU - 1.001);
  const double fv = std::min(std::max((std::log(pitchCm) - kV0) / kDV, 0.), kNV - 1.001);
  const int iu = static_cast<int>(fu), iv = static_cast<int>(fv);
  const double wu = fu - iu, wv = fv - iv;
  const double p = (1. - wu) * ((1. - wv) * CellValue(Cell(iu, iv, massMeV), dedx) +
                                wv * CellValue(Cell(iu, iv + 1, massMeV), dedx)) +
                   wu * ((1. - wv) * CellValue(Cell(iu + 1, iv, massMeV), dedx) +
                         wv * CellValue(Cell(iu + 1, iv + 1, massMeV), dedx));
  return std::log(std::max(p, kMinPdf));
}

bool ScanTLE(const std::vector<double>& rr, const std::vector<double>& dedx, const std::vector<double>& pitch,
             double trackLengthCm, double massMeV, double LmaxCm, double stepCm, double fineStepCm,
             double lowPRefineGeV, ScanResult& out) {
  out = ScanResult();
  const size_t n = rr.size();
  if (n == 0 || dedx.size() != n || (!pitch.empty() && pitch.size() != n) || !(stepCm > 0.)) return false;
  double rrMin = std::numeric_limits<double>::infinity();
  for (double r : rr)
    if (std::isfinite(r)) rrMin = std::min(rrMin, r);
  const double L0 = (std::isfinite(rrMin) && rrMin > 0.) ? -rrMin : 0.;
  if (!(LmaxCm >= L0)) return false;

  // Per-track median pitch: the angle dependence of the PDF without the hit-to-hit pitch fluctuations.
  std::vector<double> pitchUse(pitch);
  if (gModel.perHitPitch && gModel.trackMedianPitch && !pitch.empty()) {
    std::vector<double> v;
    for (double x : pitch)
      if (x > 0. && std::isfinite(x)) v.push_back(x);
    if (!v.empty()) {
      std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
      std::fill(pitchUse.begin(), pitchUse.end(), v[v.size() / 2]);
    }
  }

  auto scan = [&](double step, ScanResult& res) {
    res = ScanResult();
    for (double L = L0; L <= LmaxCm + 1e-9; L += step) {
      double ll = 0.;
      for (size_t i = 0; i < n; ++i) {
        const double range = rr[i] + L;
        if (!(range > 0.) || !std::isfinite(range)) continue;  // only the hit at the stopping point itself
        const double ke = KeFromRangeCm(range, massMeV);
        ll += LogPdf(gModel.dedxScale * dedx[i], ke, massMeV, pitchUse.empty() ? -1. : pitchUse[i]);
      }
      const double p = MomentumGeVFromRangeCm(trackLengthCm + L, massMeV);
      if (!(p > 0.) || !std::isfinite(ll)) continue;
      res.offsetCm.push_back(L);
      res.pGeV.push_back(p);
      res.logL.push_back(ll);
    }
  };

  scan(stepCm, out);
  if (out.pGeV.size() < 2) return false;
  if (fineStepCm > 0. && fineStepCm + 1e-12 < stepCm && lowPRefineGeV > 0. &&
      *std::min_element(out.pGeV.begin(), out.pGeV.end()) < lowPRefineGeV) {
    ScanResult fine;
    scan(fineStepCm, fine);
    if (fine.pGeV.size() >= 2) out = fine;
  }
  out.best = static_cast<int>(std::max_element(out.logL.begin(), out.logL.end()) - out.logL.begin());
  return true;
}

}  // namespace pdDedxPdf
