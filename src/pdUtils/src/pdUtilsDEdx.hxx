#ifndef pdUtilsDEdx_h
#define pdUtilsDEdx_h

#include "pdDataClasses.hxx"
#include <utility>

class TGraph;
class TMultiGraph;

namespace pdAnaUtils {

  struct DEdxFreeRangeFitResult {
    Float_t logLikelihood = -999.f;
    Float_t bestOffsetCm = -999.f;
    Float_t momentumGeV = -999.f;
  };

  Float_t GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG,
                            double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  Float_t GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                                   double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  Float_t dEdxLikelihood(TGraph* tg, TGraph* tg_ke, Float_t mass);

  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG,
                double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR,
                       double landauTruncMinRRCm = 0.,
                       double landauTailHitDropFraction = 0.);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRangeFit(AnaParticlePD* part, Int_t PDG, double Lmax = 500.,
                                                       double step = 0.5, double landauTruncMinRRCm = 0.,
                     double landauTailHitDropFraction = 0.);
  DEdxFreeRangeFitResult GetdEdxLikelihoodFreeRange_UpToRR_Fit(AnaParticlePD* part, Int_t PDG,
                       const double maxRR, double Lmax = 500., double step = 0.5,
                       double landauTruncMinRRCm = 0.,
                       double landauTailHitDropFraction = 0.);

    TMultiGraph* MakePionFreeRangeDedxVsRRMultiGraph(AnaParticlePD* part, double Lmax = 500., double step = 0.5,
                  double landauTruncMinRRCm = 0., double landauTailHitDropFraction = 0.,
                  const char* xAxisTitle = nullptr);

  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke, Float_t mass);
  DEdxFreeRangeFitResult dEdxLikelihoodFreeRangeFit(TGraph* tg, TGraph* tg_ke, Float_t mass, double L0,
                                                    double Lmax, double step, double measuredTrackLengthCm,
                                                    bool computeMomentum);

  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
  double KineticEnergyMeVFromResidualRangeCm(TGraph* tg_ke, double range_cm);
  double GetWmax(double KE, double mass);
  double GetLandauXi(double KE, double dx, double mass);
  double dEdxPDF(double *x, double *par);

}

#endif
