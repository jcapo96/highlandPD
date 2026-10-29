#ifndef pionMomentumDump_h
#define pionMomentumDump_h

// Per-track dump for pion momentum estimator studies (one TTree entry per reco track).
// Kinds: 0 = beam (main) track, 1 = reco daughter of the beam track that is a true charged pion (MC),
//        2 = reco track matched to a true charged pion whose true parent is a K0S (MC).
// Stores collection-plane hits, reco trajectory points, truth and the baseline TLE / MCS / CSDA estimates,
// so that estimators can be prototyped offline without rerunning the framework.

#include "pdDataClasses.hxx"
#include "pdMomReconstruction.hxx"
#include <TTree.h>
#include <vector>

class OutputManager;
class AnaEventB;

namespace pionMomentumDump {

struct DumpConfig {
  pdMomReconstruction::PionTLEFitConfig tle;
  pdMomReconstruction::PionMCSConfig mcs;
  pdMomReconstruction::PionMomentumRecoConfig reco;
  int maxTrajPoints = 2000;
  bool tleOnBeamTrack = false;  // the free-range scan is slow on long beam tracks and has no sensitivity there
  bool storeTleCurve = true;    // store the TLE lnL(p) curve (doubles the TLE time)
};

class TrackDump {
public:
  TrackDump() = default;
  ~TrackDump() = default;

  /// Fill all tracks of interest in this event. anaEntry is the entry index of the event in the ana tree.
  void FillEvent(OutputManager& output, AnaEventB& event, AnaParticlePD* mainTrack, AnaBeamPD* beam, Long64_t anaEntry,
                 const DumpConfig& cfg);
  void Write(OutputManager& output);

private:
  void Book(OutputManager& output);
  void FillTrack(AnaParticlePD* track, int kind, AnaBeamPD* beam, const DumpConfig& cfg);

  TTree* _tree = nullptr;

  Long64_t _anaentry = -1;
  Int_t _run = -999, _subrun = -999, _evt = -999, _ismc = 0, _kind = -1;
  Float_t _beammom = -999.f, _beammomtpc = -999.f;
  Int_t _beampdgpi = 0;
  // truth
  Int_t _truepdg = -999, _trueid = -999, _trueparentpdg = -999, _trueparentid = -999, _trueprocend = -999;
  Float_t _truep0 = -999.f, _truepend = -999.f, _truepintpc = -999.f, _truepatstart = -999.f, _truelength = -999.f;
  Float_t _truestart[3] = {-999.f, -999.f, -999.f}, _trueend[3] = {-999.f, -999.f, -999.f};
  Float_t _truedir[3] = {-999.f, -999.f, -999.f};
  // reco
  Float_t _length = -999.f;
  Float_t _start[3] = {-999.f, -999.f, -999.f}, _end[3] = {-999.f, -999.f, -999.f};
  Float_t _dirstart[3] = {-999.f, -999.f, -999.f};
  Int_t _ndau = 0, _ndesc = 0, _ndeschits = 0, _isbeam = 0;
  Float_t _descevis = 0.f;  // visible energy of all reco descendants, collection plane, sum dE/dx*pitch [MeV]
  Float_t _chi2proton = -999.f, _chi2ndf = -999.f;
  std::vector<float> _hrr, _hdedx, _hpitch, _hx, _hy, _hz, _htick;
  std::vector<int> _hwire, _htpc;
  std::vector<float> _tx, _ty, _tz;
  // baseline estimators [GeV/c]
  Float_t _ptle = -999.f, _pmcs = -999.f, _pcsda = -999.f;
  Float_t _prec = -999.f;  // EstimatePionMomentumRecommended
  Int_t _precmeasured = 0;
  std::vector<float> _tlep, _tlell;
};

}  // namespace pionMomentumDump

#endif
