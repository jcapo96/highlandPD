// Class:       PointIdAlg
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),      from DUNE,   CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),            from DUNE,   FNAL/NCBJ, since May 2016
//              P.Plonski,                                 from DUNE,   WUT,       since May 2016
//              D.Smith,                                   from LArIAT, BU, 2017: real data dump

#include "PointIdAlg.hxx"

#ifdef  CompileTF
#include "tensorflow/core/public/session.h"  // anselmo
#endif
/*
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireIDError
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larsim/Simulation/LArG4Parameters.h"

#include "CLHEP/Random/RandGauss.h"
*/
#include <sys/stat.h>


std::string wspaces2(size_t n){
  std::string ws="";
  for (size_t i=0;i<n;i++)
    ws +=" ";
  return ws;
}



// ------------------------------------------------------
// -------------------ModelInterface---------------------
// ------------------------------------------------------

std::vector<std::vector<float>>
nnet::ModelInterface::Run(std::vector<std::vector<std::vector<float>>> const& inps, int samples)
{
  if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty())
    return std::vector<std::vector<float>>();

  if ((samples == -1) || (samples > (int)inps.size())) { samples = inps.size(); }

  std::vector<std::vector<float>> results;
  for (int i = 0; i < samples; ++i) {
    results.push_back(Run(inps[i]));
  }
  return results;
}

std::string
nnet::ModelInterface::findFile(const char* fileName) const
{
  std::string fname_out;
  std::string path = "";
  /*
  cet::search_path sp("FW_SEARCH_PATH");
  if (!sp.find_file(fileName, fname_out)) {
    struct stat buffer;
    if (stat(fileName, &buffer) == 0) { fname_out = fileName; }
    else {
      throw art::Exception(art::errors::NotFound) << "Could not find the model file " << fileName;
    }
  }
  */
  fname_out =  path + fileName;
  std::cout << fname_out << std::endl;
  return fname_out;
}

/*
// ------------------------------------------------------
// ----------------KerasModelInterface-------------------
// ------------------------------------------------------

nnet::KerasModelInterface::KerasModelInterface(const char* modelFileName)
  : m(nnet::ModelInterface::findFile(modelFileName).c_str())
{
  mf::LogInfo("KerasModelInterface") << "Keras model loaded.";
}
// ------------------------------------------------------

std::vector<float>
nnet::KerasModelInterface::Run(std::vector<std::vector<float>> const& inp2d)
{
  std::vector<std::vector<std::vector<float>>> inp3d;
  inp3d.push_back(inp2d); // lots of copy, should add 2D to keras...

  keras::DataChunk* sample = new keras::DataChunk2D();
  sample->set_data(inp3d); // and more copy...
  std::vector<float> out = m.compute_output(sample);
  delete sample;
  return out;
}
*/
// ------------------------------------------------------
// -----------------TfModelInterface---------------------
// ------------------------------------------------------

nnet::TfModelInterface::TfModelInterface(const char* modelFileName)
{

#ifdef  CompileTF
  g = tf::Graph::create(nnet::ModelInterface::findFile(modelFileName).c_str(),
                        {"cnn_output", "_netout"});
  //if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
  if (!g) {std::cout << "TF model failed." << std::endl; exit(1);}
#endif
  std::cout << "TF model loaded." << std::endl;
  
}
// ------------------------------------------------------

std::vector<std::vector<float>>
nnet::TfModelInterface::Run(std::vector<std::vector<std::vector<float>>> const& inps, int samples)
{
  if (debug_level>=2) std::cout << wspaces2(4) << "PointIdAlg.cxx:Run. inps.size() = " << inps.size() << std::endl;   
  if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty())
    return std::vector<std::vector<float>>();

  if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

  long long int rows = inps.front().size(), cols = inps.front().front().size();

  if (debug_level>=3) std::cout << wspaces2(6) << "PointIdAlg.cxx:Run. rows, cols = " << rows << " " << cols << std::endl;   
  
  if (debug_level>=3) std::cout << wspaces2(6) << "PointIdAlg.cxx:Run. READY TO CALL TF " << std::endl;   

  size_t count_gt0=0;

#ifdef  CompileTF
  tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({samples, rows, cols, 1}));   // anselmo
  auto input_map = _x.tensor<float, 4>();
#endif
  for (long long int s = 0; s < samples; ++s) {
    const auto& sample = inps[s];
    for (long long int r = 0; r < rows; ++r) {
      const auto& row = sample[r];
      for (long long int c = 0; c < cols; ++c) {
#ifdef  CompileTF
        input_map(s, r, c, 0) = row[c];  // anselmo
#endif
        if (row[c]>0) count_gt0++;
        if (debug_level>=5) std::cout << wspaces2(8) << "PointIdAlg.cxx:Run. s, r, c, row[c] = " << s << " " << r << " " << c << " " << row[c] << std::endl;
        else if (debug_level>=4 && row[c]>0 ) std::cout << wspaces2(6) << "PointIdAlg.cxx:Run. s, r, c, row[c] = " << s << " " << r << " " << c << " " << row[c] << std::endl;   
      }
    }
  }

  if (debug_level>=3) std::cout << wspaces2(6) << "PointIdAlg.cxx:Run. #tensor elements > 0 = " << count_gt0 << std::endl;   

#ifdef  CompileTF 
  return g->run(_x);
#else
  return std::vector<std::vector<float>>();
#endif
}
// ------------------------------------------------------

std::vector<float>
nnet::TfModelInterface::Run(std::vector<std::vector<float>> const& inp2d)
{
  long long int rows = inp2d.size(), cols = inp2d.front().size();

  if (debug_level>=2) std::cout << wspaces2(4) << "PointIdAlg.cxx:Run. rows, cols = " << rows << " " << cols << std::endl;   
  /*
  tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, rows, cols, 1}));
  auto input_map = _x.tensor<float, 4>();
  for (long long int r = 0; r < rows; ++r) {
    const auto& row = inp2d[r];
    for (long long int c = 0; c < cols; ++c) {
      input_map(0, r, c, 0) = row[c];
    }
  }

  auto out = g->run(_x);
  if (!out.empty())
    return out.front();
  else
    return std::vector<float>();
  */
  return std::vector<float>();
}

// ------------------------------------------------------
// --------------------PointIdAlg------------------------
// ------------------------------------------------------

  // anselmo: Taken from https://internal.dunescience.org/doxygen/standard__reco__dune10kt_8fcl_source.html
nnet::PointIdAlg::PointIdAlg()
  : // img::DataProviderAlg(config)
  fNNet(0)
  , fPatchSizeW(48)
  , fPatchSizeD(48)
  , fCurrentWireIdx(99999)
  , fCurrentScaledDrift(99999)
{
  fNNetModelFilePath = "cnn.pb";
  //  fNNetOutputs = "";
  
  deleteNNet();
  /*
  if ((fNNetModelFilePath.length() > 5) &&
      (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 5, 5, ".nnet") == 0)) {
    fNNet = new nnet::KerasModelInterface(fNNetModelFilePath.c_str());
  }
  */    
  if ((fNNetModelFilePath.length() > 3) &&
      (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 3, 3, ".pb") == 0)) {    
    fNNet = new nnet::TfModelInterface(fNNetModelFilePath.c_str());
  }
  else {
    std::cout  << "File name extension not supported. Only .pb files are supported" << std::endl;
  }

  if (!fNNet) { std:: cout << "Loading model from file failed." << std::endl;exit(1);}

  resizePatch();
}

// ------------------------------------------------------

nnet::PointIdAlg::~PointIdAlg()
{
  deleteNNet();
}
// ------------------------------------------------------

void
nnet::PointIdAlg::resizePatch()
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx::resizePatch" << std::endl;   
  fWireDriftPatch.resize(fPatchSizeW);
  for (auto& r : fWireDriftPatch)
    r.resize(fPatchSizeD);
}
// ------------------------------------------------------

float
nnet::PointIdAlg::predictIdValue(unsigned int wire, float drift, size_t outIdx) const
{
  if (debug_level>=1) std::cout << wspaces2(2) << "PointIdAlg.cxx:predictIdValue" << std::endl;   
  float result = 0.;

  if (!bufferPatch(wire, drift)) {
    //    mf::LogError("PointIdAlg") << "Patch buffering failed.";
    std::cout << "Patch buffering failed." << std::endl;
    return result;
  }

  if (fNNet) {
    auto out = fNNet->Run(fWireDriftPatch);
    if (!out.empty()) { result = out[outIdx]; }
    else {
      //      mf::LogError("PointIdAlg") << "Problem with applying model to input.";
      std::cout << "Problem with applying model to input." << std::endl;
    }
  }

  return result;
}
// ------------------------------------------------------

std::vector<float>
nnet::PointIdAlg::predictIdVector(unsigned int wire, float drift) const
{
  if (debug_level>=1) std::cout << wspaces2(2) << "PointIdAlg.cxx:predictIdVector" << std::endl;   
  std::vector<float> result;
  if (!bufferPatch(wire, drift)) {
    //    mf::LogError("PointIdAlg") << "Patch buffering failed.";
    std::cout << "Patch buffering failed." << std::endl;
    return result;
  }

  if (fNNet) {
    result = fNNet->Run(fWireDriftPatch);
    if (result.empty()) //{ mf::LogError("PointIdAlg") << "Problem with applying model to input."; }
      { std::cout << "Problem with applying model to input." << std::endl; }
  }

  return result;
}
// ------------------------------------------------------

std::vector<std::vector<float>>
nnet::PointIdAlg::predictIdVectors(std::vector<std::pair<unsigned int, float>> points) const
{
  if (debug_level>=1) std::cout << wspaces2(2) << "PointIdAlg.cxx:predictIdVectors" << std::endl;   
  if (points.empty() || !fNNet) { return std::vector<std::vector<float>>(); }
  std::vector<std::vector<std::vector<float>>> inps(
                                                    points.size(), std::vector<std::vector<float>>(fPatchSizeW, std::vector<float>(fPatchSizeD)));
  for (size_t i = 0; i < points.size(); ++i) {
    unsigned int wire = points[i].first;
    float drift = points[i].second;
    if (!bufferPatch(wire, drift, inps[i])) {
      //      throw cet::exception("PointIdAlg") << "Patch buffering failed" << std::endl;
      std::cout << "Patch buffering failed" << std::endl;
    }    
  }
  
  return fNNet->Run(inps);
}
// ------------------------------------------------------

bool
nnet::PointIdAlg::isSamePatch(unsigned int wire1,
                              float drift1,
                              unsigned int wire2,
                              float drift2) const
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx:isSamePatch" << std::endl;   
  if (fDownscaleFullView) {
    size_t sd1 = (size_t)(drift1 / fDriftWindow);
    size_t sd2 = (size_t)(drift2 / fDriftWindow);
    if ((wire1 == wire2) && (sd1 == sd2)) return true; // the same position
  }
  else {
    if ((wire1 == wire2) && ((size_t)drift1 == (size_t)drift2)) return true; // the same position
  }

  return false; // not the same position
}

bool
nnet::PointIdAlg::isCurrentPatch(unsigned int wire, float drift) const
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx:isCurrentPatch" << std::endl;   
  if (fDownscaleFullView) {
    size_t sd = (size_t)(drift / fDriftWindow);
    if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == sd))
      return true; // still within the current position
  }
  else {
    if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == drift))
      return true; // still within the current position
  }

  return false; // not a current position
}
// ------------------------------------------------------

std::vector<float>
nnet::PointIdAlg::flattenData2D(std::vector<std::vector<float>> const& patch)
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx: flattenData2D " << std::endl;   
  std::vector<float> flat;
  if (patch.empty() || patch.front().empty()) {
    //    mf::LogError("DataProviderAlg") << "Patch is empty.";
    std::cout << "Patch is empty." << std::endl;
    return flat;
  }

  flat.resize(patch.size() * patch.front().size());

  for (size_t w = 0, i = 0; w < patch.size(); ++w) {
    auto const& wire = patch[w];
    for (size_t d = 0; d < wire.size(); ++d, ++i) {
      flat[i] = wire[d];
    }
  }

  return flat;
}
// ------------------------------------------------------

bool
nnet::PointIdAlg::isInsideFiducialRegion(unsigned int wire, float drift) const
{
  if (debug_level>=2) std::cout << wspaces2(4) << "PointIdAlg.cxx:isInsideFiducialRegion" << std::endl;   
  size_t marginW = fPatchSizeW / 8; // fPatchSizeX/2 will make patch always completely filled
  size_t marginD = fPatchSizeD / 8;

  size_t scaledDrift = (size_t)(drift / fDriftWindow);
  if ((wire >= marginW) && (wire < fAlgView.fNWires - marginW) && (scaledDrift >= marginD) &&
      (scaledDrift < fAlgView.fNScaledDrifts - marginD))
    return true;
  else
    return false;
}
// ------------------------------------------------------

// ------------------------------------------------------
// ------------------TrainingDataAlg---------------------
// ------------------------------------------------------
/*
nnet::TrainingDataAlg::TrainingDataAlg(const Config& config)
  : img::DataProviderAlg(config)
  , fEdepTot(0)
  , fWireProducerLabel(config.WireLabel())
  , fHitProducerLabel(config.HitLabel())
  , fTrackModuleLabel(config.TrackLabel())
  , fSimulationProducerLabel(config.SimulationLabel())
  , fSaveVtxFlags(config.SaveVtxFlags())
  , fAdcDelay(config.AdcDelayTicks())
  , fEventsPerBin(100, 0)
{
  fSaveSimInfo = !fSimulationProducerLabel.label().empty();
}
*/
// ------------------------------------------------------

nnet::TrainingDataAlg::~TrainingDataAlg() = default;
// ------------------------------------------------------

img::DataProviderAlgView
nnet::TrainingDataAlg::resizeView(detinfo::DetectorClocksData const& clock_data,
                                  detinfo::DetectorPropertiesData const& det_prop,
                                  size_t wires,
                                  size_t drifts)
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx:resizeView" << std::endl;   
  auto view = img::DataProviderAlg::resizeView(clock_data, det_prop, wires, drifts);

  fWireDriftEdep.resize(wires);
  for (auto& w : fWireDriftEdep) {
    w.resize(view.fNCachedDrifts);
    std::fill(w.begin(), w.end(), 0.0F);
  }

  fWireDriftPdg.resize(wires);
  for (auto& w : fWireDriftPdg) {
    w.resize(view.fNCachedDrifts);
    std::fill(w.begin(), w.end(), 0);
  }
  return view;
}
// ------------------------------------------------------

bool
nnet::TrainingDataAlg::setWireEdepsAndLabels(std::vector<float> const& edeps,
                                             std::vector<int> const& pdgs,
                                             size_t wireIdx)
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx:setWireEdepsAndLabels" << std::endl;     
  if ((wireIdx >= fWireDriftEdep.size()) || (edeps.size() != pdgs.size())) { return false; }

  size_t dstep = 1;
  if (fDownscaleFullView) { dstep = fDriftWindow; }

  if (edeps.size() / dstep > fAlgView.fNCachedDrifts) { return false; }

  auto& wEdep = fWireDriftEdep[wireIdx];
  auto& wPdg = fWireDriftPdg[wireIdx];

  for (size_t i = 0; i < fAlgView.fNCachedDrifts; ++i) {
    size_t i0 = i * dstep;
    size_t i1 = (i + 1) * dstep;

    int best_pdg = pdgs[i0] & nnet::TrainingDataAlg::kPdgMask;
    int vtx_flags = pdgs[i0] & nnet::TrainingDataAlg::kVtxMask;
    int type_flags = pdgs[i0] & nnet::TrainingDataAlg::kTypeMask;
    float max_edep = edeps[i0];
    for (size_t k = i0 + 1; k < i1; ++k) {
      float ek = edeps[k];
      if (ek > max_edep) {
        max_edep = ek;
        best_pdg = pdgs[k] & nnet::TrainingDataAlg::kPdgMask; // remember best matching pdg
      }
      type_flags |= pdgs[k] & nnet::TrainingDataAlg::kTypeMask; // accumulate track type flags
      vtx_flags |= pdgs[k] & nnet::TrainingDataAlg::kVtxMask;   // accumulate all vtx flags
    }

    wEdep[i] = max_edep;

    best_pdg |= type_flags;
    if (fSaveVtxFlags) best_pdg |= vtx_flags;
    wPdg[i] = best_pdg;
  }

  return true;
}

// ------------------------------------------------------

nnet::TrainingDataAlg::WireDrift
nnet::TrainingDataAlg::getProjection(detinfo::DetectorClocksData const& clockData,
                                     detinfo::DetectorPropertiesData const& detProp,
                                     const TLorentzVector& tvec,
                                     unsigned int plane) const
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx:getProjection" << std::endl;     
  nnet::TrainingDataAlg::WireDrift wd;
  wd.Wire = 0;
  wd.Drift = 0;
  wd.TPC = -1;
  wd.Cryo = -1;
  /*
  try {
    double vtx[3] = {tvec.X(), tvec.Y(), tvec.Z()};
    if (fGeometry->FindTPCAtPosition(vtx).isValid) {
      geo::TPCID tpcid = fGeometry->FindTPCAtPosition(vtx);
      unsigned int tpc = tpcid.TPC, cryo = tpcid.Cryostat;

      // correct for the time offset
      float dx = tvec.T() * 1.e-3 * detProp.DriftVelocity();
      int driftDir = fGeometry->TPC(tpcid).DetectDriftDirection();
      if (driftDir == 1) { dx *= -1; }
      else if (driftDir != -1) {
        throw cet::exception("nnet::TrainingDataAlg") << "drift direction is not X." << std::endl;
      }
      vtx[0] = tvec.X() + dx;

      wd.Wire = fGeometry->NearestWire(vtx, plane, tpc, cryo);
      wd.Drift = detProp.ConvertXToTicks(vtx[0], plane, tpc, cryo);
      wd.TPC = tpc;
      wd.Cryo = cryo;
    }
  }
  catch (const geo::InvalidWireIDError& e) {
    mf::LogWarning("TrainingDataAlg")
      << "Vertex projection out of wire planes, just skipping this vertex.";
  }
  catch (...) {
    mf::LogWarning("TrainingDataAlg") << "Vertex projection out of wire planes, skip MC vertex.";
  }
  */
  return wd;
}
// ------------------------------------------------------
/*
bool
nnet::TrainingDataAlg::isElectronEnd(
  const simb::MCParticle& particle,
  const std::unordered_map<int, const simb::MCParticle*>& particleMap) const
{
  const float minElectronLength2 = 2.5 * 2.5;
  const float maxDeltaLength2 = 0.15 * 0.15;

  int pdg = abs(particle.PdgCode());
  if (pdg != 11) return false; // should be applied only to electrons

  size_t nSec = particle.NumberDaughters();
  for (size_t d = 0; d < nSec; ++d) {
    auto d_search = particleMap.find(particle.Daughter(d));
    if (d_search != particleMap.end()) {
      auto const& daughter = *((*d_search).second);
      int d_pdg = abs(daughter.PdgCode());
      if (d_pdg != 22) { return false; } // not the end of the shower
    }
  }

  float trkLength2 = 0;
  auto const* p = &particle;
  bool branching = false;
  while (!branching) {
    trkLength2 += particleRange2(*p);
    auto m_search = particleMap.find(p->Mother());
    if (m_search != particleMap.end()) {
      p = (*m_search).second;
      int m_pdg = abs(p->PdgCode());
      if (m_pdg == 11) {
        nSec = p->NumberDaughters();
        size_t ne = 0;
        for (size_t d = 0; d < nSec; ++d) {
          auto d_search = particleMap.find(p->Daughter(d));
          if (d_search != particleMap.end()) {
            auto const& daughter = *((*d_search).second);
            int d_pdg = abs(daughter.PdgCode());
            if (d_pdg == 11) {
              if (particleRange2(daughter) > maxDeltaLength2) { ne++; }
            }
          }
        }
        if (ne > 1) { branching = true; }
      }
      else
        break;
    }
    else
      break;
  }

  return (trkLength2 > minElectronLength2);
}

bool
nnet::TrainingDataAlg::isMuonDecaying(
  const simb::MCParticle& particle,
  const std::unordered_map<int, const simb::MCParticle*>& particleMap) const
{
  bool hasElectron = false, hasNuMu = false, hasNuE = false;

  int pdg = abs(particle.PdgCode());
  if ((pdg == 13) && (particle.EndProcess() == "FastScintillation")) // potential muon decay at rest
  {
    unsigned int nSec = particle.NumberDaughters();
    for (size_t d = 0; d < nSec; ++d) {
      auto d_search = particleMap.find(particle.Daughter(d));
      if (d_search != particleMap.end()) {
        auto const& daughter = *((*d_search).second);
        int d_pdg = abs(daughter.PdgCode());
        if (d_pdg == 11)
          hasElectron = true;
        else if (d_pdg == 14)
          hasNuMu = true;
        else if (d_pdg == 12)
          hasNuE = true;
      }
    }
  }

  return (hasElectron && hasNuMu && hasNuE);
}

void
nnet::TrainingDataAlg::collectVtxFlags(
  std::unordered_map<size_t, std::unordered_map<int, int>>& wireToDriftToVtxFlags,
  detinfo::DetectorClocksData const& clockData,
  detinfo::DetectorPropertiesData const& detProp,
  const std::unordered_map<int, const simb::MCParticle*>& particleMap,
  unsigned int plane) const
{
  for (auto const& p : particleMap) {
    auto const& particle = *p.second;

    double ekStart = 1000. * (particle.E() - particle.Mass());
    double ekEnd = 1000. * (particle.EndE() - particle.Mass());
    int pdg = abs(particle.PdgCode());
    int flagsStart = nnet::TrainingDataAlg::kNone;
    int flagsEnd = nnet::TrainingDataAlg::kNone;

    switch (pdg) {
    case 22:                                                     // gamma
      if ((particle.EndProcess() == "conv") && (ekStart > 40.0)) // conversion, gamma > 40MeV
      {
        flagsEnd = nnet::TrainingDataAlg::kConv;
      }
      break;

    case 11: // e+/-
      if (isElectronEnd(particle, particleMap)) { flagsEnd = nnet::TrainingDataAlg::kElectronEnd; }
      break;

    case 13: // mu+/-
      if (isMuonDecaying(particle, particleMap)) { flagsEnd = nnet::TrainingDataAlg::kDecay; }
      break;

    case 111: // pi0
      flagsStart = nnet::TrainingDataAlg::kPi0;
      break;

    case 321:  // K+/-
    case 211:  // pi+/-
    case 2212: // proton
      if (ekStart > 50.0) {
        if (particle.Mother() != 0) {
          auto search = particleMap.find(particle.Mother());
          if (search != particleMap.end()) {
            auto const& mother = *((*search).second);
            int m_pdg = abs(mother.PdgCode());
            unsigned int nSec = mother.NumberDaughters();
            unsigned int nVisible = 0;
            if (nSec > 1) {
              for (size_t d = 0; d < nSec; ++d) {
                auto d_search = particleMap.find(mother.Daughter(d));
                if (d_search != particleMap.end()) {
                  auto const& daughter = *((*d_search).second);
                  int d_pdg = abs(daughter.PdgCode());
                  if (((d_pdg == 2212) || (d_pdg == 211) || (d_pdg == 321)) &&
                      (1000. * (daughter.E() - daughter.Mass()) > 50.0)) {
                    ++nVisible;
                  }
                }
              }
            }
            // hadron with Ek > 50MeV (so well visible) and
            // produced by another hadron (but not neutron, so not single track from nothing) or
            // at least secondary hadrons with Ek > 50MeV (so this is a good kink or V-like)
            if (((m_pdg != pdg) && (m_pdg != 2112)) || ((m_pdg != 2112) && (nVisible > 0)) ||
                ((m_pdg == 2112) && (nVisible > 1))) {
              flagsStart = nnet::TrainingDataAlg::kHadr;
            }
          }
        }

        if (particle.EndProcess() == "FastScintillation") // potential decay at rest
        {
          unsigned int nSec = particle.NumberDaughters();
          for (size_t d = 0; d < nSec; ++d) {
            auto d_search = particleMap.find(particle.Daughter(d));
            if (d_search != particleMap.end()) {
              auto const& daughter = *((*d_search).second);
              int d_pdg = abs(daughter.PdgCode());
              if ((pdg == 321) && (d_pdg == 13)) {
                flagsEnd = nnet::TrainingDataAlg::kDecay;
                break;
              }
              if ((pdg == 211) && (d_pdg == 13)) {
                flagsEnd = nnet::TrainingDataAlg::kDecay;
                break;
              }
            }
          }
        }

        if ((particle.EndProcess() == "Decay") && (ekEnd > 200.0)) // decay in flight
        {
          unsigned int nSec = particle.NumberDaughters();
          for (size_t d = 0; d < nSec; ++d) {
            auto d_search = particleMap.find(particle.Daughter(d));
            if (d_search != particleMap.end()) {
              auto const& daughter = *((*d_search).second);
              int d_pdg = abs(daughter.PdgCode());
              if ((pdg == 321) && (d_pdg == 13)) {
                flagsEnd = nnet::TrainingDataAlg::kHadr;
                break;
              }
              if ((pdg == 211) && (d_pdg == 13)) {
                flagsEnd = nnet::TrainingDataAlg::kHadr;
                break;
              }
            }
          }
        }
      }
      break;

    default: continue;
    }

    if (particle.Process() == "primary") { flagsStart |= nnet::TrainingDataAlg::kNuPri; }

    if (flagsStart != nnet::TrainingDataAlg::kNone) {
      auto wd = getProjection(clockData, detProp, particle.Position(), plane);

      if ((wd.TPC == TPC()) && (wd.Cryo == Cryo())) {
        wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsStart;
      }
    }
    if (flagsEnd != nnet::TrainingDataAlg::kNone) {
      auto wd = getProjection(clockData, detProp, particle.EndPosition(), plane);
      if ((wd.TPC == TPC()) && (wd.Cryo == Cryo())) {
        wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= flagsEnd;
      }
    }

    //TY: check elastic/inelastic scattering
    if (pdg == 321 || pdg == 211 || pdg == 2212) {
      simb::MCTrajectory truetraj = particle.Trajectory();
      auto thisTrajectoryProcessMap1 = truetraj.TrajectoryProcesses();
      if (thisTrajectoryProcessMap1.size()) {
        for (auto const& couple1 : thisTrajectoryProcessMap1) {
          if ((truetraj.KeyToProcess(couple1.second)).find("Elastic") != std::string::npos) {
            auto wd = getProjection(clockData, detProp, truetraj.at(couple1.first).first, plane);
            if ((wd.TPC == TPC()) && (wd.Cryo == Cryo())) {
              wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= nnet::TrainingDataAlg::kElastic;
            }
          }
          if ((truetraj.KeyToProcess(couple1.second)).find("Inelastic") != std::string::npos) {
            auto wd = getProjection(clockData, detProp, truetraj.at(couple1.first).first, plane);
            if ((wd.TPC == TPC()) && (wd.Cryo == Cryo())) {
              wireToDriftToVtxFlags[wd.Wire][wd.Drift] |= nnet::TrainingDataAlg::kInelastic;
            }
          }
        }
      }
    }
  }
}
*/


// ------------------------------------------------------

bool
nnet::TrainingDataAlg::setDataEventData(const AnaEvent& event,
                                        detinfo::DetectorClocksData const& clockData,
                                        detinfo::DetectorPropertiesData const& detProp,
                                        unsigned int plane,
                                        unsigned int tpc,
                                        unsigned int cryo)
{

  if (debug_level>=0) std::cout << "PointIdAlg.cxx: setDataEventData " << std::endl;   
  std::vector<my_recob::Wire> wireHandle;// art::Handle<std::vector<my_recob::Wire>> wireHandle; 
  std::vector<my_recob::Wire*> Wirelist; // std::vector<art::Ptr<my_recob::Wire>> Wirelist;      

  //  if (event.getByLabel(fWireProducerLabel, wireHandle)) art::fill_ptr_vector(Wirelist, wireHandle);

  if (!setWireDriftData(clockData, detProp, wireHandle, plane, tpc, cryo)) {
    //    mf::LogError("TrainingDataAlg") << "Wire data not set.";
    std::cout << "Wire data not set." << std::endl;;
    return false;
  }

  // Hit info
  std::vector<AnaHitPD*> HitHandle;  // art::Handle<std::vector<my_recob::Hit>> HitHandle;    
  std::vector<AnaHitPD*> Hitlist;    // std::vector<art::Ptr<my_recob::Hit>> Hitlist;         

  //  if (event.getByLabel(fHitProducerLabel, HitHandle)) art::fill_ptr_vector(Hitlist, HitHandle);

  // Track info
  std::vector<my_recob::Track> TrackHandle;// art::Handle<std::vector<my_recob::Track>> TrackHandle;
  std::vector<my_recob::Track> Tracklist;  // std::vector<art::Ptr<my_recob::Track>> Tracklist;  

  //  if (event.getByLabel(fTrackModuleLabel, TrackHandle))
  //    art::fill_ptr_vector(Tracklist, TrackHandle);

  //  art::FindManyP<my_recob::Track> ass_trk_hits(HitHandle, event, fTrackModuleLabel);
  std::vector<std::vector<my_recob::Track*> >ass_trk_hits;
  
  // Loop over wires (sorry about hard coded value) to fill in 1) pdg and 2) charge depo
  for (size_t widx = 0; widx < 240; ++widx) {

    std::vector<float> labels_deposit(fAlgView.fNDrifts, 0); // full-drift-length buffers
    std::vector<int> labels_pdg(fAlgView.fNDrifts, 0);

    // First, the charge depo
    for (size_t subwidx = 0; subwidx < Wirelist.size(); ++subwidx) {
      if (widx + 240 == Wirelist[subwidx]->Channel()) {
        labels_deposit = Wirelist[subwidx]->Signal();
        break;
      }
    }

    // Second, the pdg code
    // This code finds the angle of the track and records
    //  events based on its angle to try to get an isometric sample
    //  instead of just a bunch of straight tracks

    // Meta code:
    // For each hit:
    //  find farthest hit from point
    //  then find farthest hit from THAT one
    //  should be start and end of track, then just use trig

    for (size_t iHit = 0; iHit < Hitlist.size(); ++iHit) {

      if (Hitlist[iHit]->Channel != widx + 240) { continue; }
      if (Hitlist[iHit]->View != 1) { continue; }

      // Make sure there is a track association
      if (ass_trk_hits.at(iHit).size() == 0) { continue; }

      // Not sure about this
      // Cutting on length to not just get a bunch of shower stubs
      // Might add a lot of bias though
      if (ass_trk_hits.at(iHit)[0]->Length() < 5) { continue; }

      // Search for farest hit from this one
      int far_index = 0;
      double far_dist = 0;

      for (size_t jHit = 0; jHit < Hitlist.size(); ++jHit) {
        if (jHit == iHit) { continue; }
        if (Hitlist[jHit]->View != 1) { continue; }

        if (ass_trk_hits.at(jHit).size() == 0) { continue; }
        if (ass_trk_hits.at(jHit)[0]->ID() != ass_trk_hits.at(iHit)[0]->ID()) { continue; }

        double dist = sqrt((Hitlist[iHit]->Channel - Hitlist[jHit]->Channel) *
                             (Hitlist[iHit]->Channel - Hitlist[jHit]->Channel) +
                           (Hitlist[iHit]->PeakTime - Hitlist[jHit]->PeakTime) *
                             (Hitlist[iHit]->PeakTime - Hitlist[jHit]->PeakTime));

        if (far_dist < dist) {
          far_dist = dist;
          far_index = jHit;
        }
      }

      // Search for the other end of the track
      int other_end = 0;
      int other_dist = 0;

      for (size_t jHit = 0; jHit < Hitlist.size(); ++jHit) {
        if (jHit == iHit or int(jHit) == far_index) { continue; }
        if (Hitlist[jHit]->View != 1) { continue; }

        if (ass_trk_hits.at(jHit).size() == 0) { continue; }
        if (ass_trk_hits.at(jHit)[0]->ID() != ass_trk_hits.at(iHit)[0]->ID()) { continue; }

        double dist = sqrt((Hitlist[far_index]->Channel - Hitlist[jHit]->Channel) *
                             (Hitlist[far_index]->Channel - Hitlist[jHit]->Channel) +
                           (Hitlist[far_index]->PeakTime - Hitlist[jHit]->PeakTime) *
                             (Hitlist[far_index]->PeakTime - Hitlist[jHit]->PeakTime));

        if (other_dist < dist) {
          other_dist = dist;
          other_end = jHit;
        }
      }

      // We have the end points now
      double del_wire = double(Hitlist[other_end]->Channel - Hitlist[far_index]->Channel);
      double del_time = double(Hitlist[other_end]->PeakTime - Hitlist[far_index]->PeakTime);
      double hypo = sqrt(del_wire * del_wire + del_time * del_time);

      if (hypo == 0) { continue; } // Should never happen, but doing it anyway

      double cosser = TMath::Abs(del_wire / hypo);
      double norm_ang = TMath::ACos(cosser) * 2 / TMath::Pi();

      // Using fEventsPerBin to keep track of number of hits per angle (normalized to 0 to 1)

      int binner = int(norm_ang * fEventsPerBin.size());
      if (binner >= (int)fEventsPerBin.size()) {
        binner = fEventsPerBin.size() - 1;
      } // Dealing with rounding errors

      // So we should get a total of 5000 * 100 = 50,000 if we use the whole set
      if (fEventsPerBin[binner] > 5000) { continue; }
      fEventsPerBin[binner]++;

      // If survives everything, saves the pdg
      labels_pdg[Hitlist[iHit]->PeakTime] = 211; // Same as pion for now
    }

    setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);

  } // for each Wire

  return true;
}



// ------------------------------------------------------

/*
bool
nnet::TrainingDataAlg::setEventData(const AnaEvent& event,
                                    detinfo::DetectorClocksData const& clockData,
                                    detinfo::DetectorPropertiesData const& detProp,
                                    unsigned int plane,
                                    unsigned int tpc,
                                    unsigned int cryo)
{
  art::ValidHandle<std::vector<my_recob::Wire>> wireHandle =
    event.getValidHandle<std::vector<my_recob::Wire>>(fWireProducerLabel);

  if (!setWireDriftData(clockData, detProp, *wireHandle, plane, tpc, cryo)) {
    mf::LogError("TrainingDataAlg") << "Wire data not set.";
    return false;
  }

  if (!fSaveSimInfo || event.isRealData()) {
    mf::LogInfo("TrainingDataAlg") << "Skip MC simulation info.";
    return true;
  }

  art::ServiceHandle<sim::LArG4Parameters const> larParameters;
  double electronsToGeV = 1. / larParameters->GeVToElectrons();

  auto particleHandle =
    event.getValidHandle<std::vector<simb::MCParticle>>(fSimulationProducerLabel);

  auto simChannelHandle =
    event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

  std::unordered_map<int, const simb::MCParticle*> particleMap;
  for (auto const& particle : *particleHandle) {
    particleMap[particle.TrackId()] = &particle;
  }

  std::unordered_map<size_t, std::unordered_map<int, int>> wireToDriftToVtxFlags;
  if (fSaveVtxFlags) collectVtxFlags(wireToDriftToVtxFlags, clockData, detProp, particleMap, plane);

  fEdepTot = 0;

  std::map<int, int> trackToPDG;
  for (size_t widx = 0; widx < fAlgView.fNWires; ++widx) {
    auto wireChannelNumber = fAlgView.fWireChannels[widx];
    if (wireChannelNumber == raw::InvalidChannelID) continue;

    std::vector<float> labels_deposit(fAlgView.fNDrifts, 0); // full-drift-length buffers,
    std::vector<int> labels_pdg(labels_deposit.size(), 0);   // both of the same size,
    int labels_size = labels_deposit.size();                 // cached as int for comparisons below

    std::map<int, std::map<int, double>> timeToTrackToCharge;
    for (auto const& channel : *simChannelHandle) {
      if (channel.Channel != wireChannelNumber) continue;

      auto const& timeSlices = channel.TDCIDEMap();
      for (auto const& timeSlice : timeSlices) {
        int time = timeSlice.first;

        auto const& energyDeposits = timeSlice.second;
        for (auto const& energyDeposit : energyDeposits) {
          int pdg = 0;
          int tid = energyDeposit.trackID;
          if (tid < 0) // negative tid means it is EM activity, and -tid is the mother
          {
            pdg = 11;
            tid = -tid;

            auto search = particleMap.find(tid);
            if (search == particleMap.end()) {
              mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
              continue;
            }
            auto const& mother = *((*search).second); // mother particle of this EM
            int mPdg = abs(mother.PdgCode());
            if ((mPdg == 13) || (mPdg == 211) || (mPdg == 2212)) {
              if (energyDeposit.numElectrons > 10)
                pdg |= nnet::TrainingDataAlg::kDelta; // tag delta ray
            }
          }
          else {
            auto search = particleMap.find(tid);
            if (search == particleMap.end()) {
              mf::LogWarning("TrainingDataAlg") << "PARTICLE NOT FOUND";
              continue;
            }
            auto const& particle = *((*search).second);
            pdg = abs(particle.PdgCode());

            if (particle.Process() == "primary") {
              if (pdg == 11) {
                pdg |= nnet::TrainingDataAlg::kPriEl; // tag primary
              }
              else if (pdg == 13) {
                pdg |= nnet::TrainingDataAlg::kPriMu; // tag primary
              }
            }

            auto msearch = particleMap.find(particle.Mother());
            if (msearch != particleMap.end()) {
              auto const& mother = *((*msearch).second);
              if (pdg == 11) // electron, check if it is Michel
              {
                if (nnet::TrainingDataAlg::isMuonDecaying(mother, particleMap)) {
                  pdg |= nnet::TrainingDataAlg::kMichel; // tag Michel
                }
              }
            }
          }

          trackToPDG[energyDeposit.trackID] = pdg;

          double energy = energyDeposit.numElectrons * electronsToGeV;
          timeToTrackToCharge[time][energyDeposit.trackID] += energy;
          fEdepTot += energy;

        } // loop over energy deposits
      }   // loop over time slices
    }     // for each SimChannel

    int type_pdg_mask = nnet::TrainingDataAlg::kTypeMask | nnet::TrainingDataAlg::kPdgMask;
    for (auto const& ttc : timeToTrackToCharge) {
      float max_deposit = 0.0;
      int max_pdg = 0;
      for (auto const& tc : ttc.second) {

        if (tc.second > max_deposit) {
          max_deposit = tc.second;
          max_pdg = trackToPDG[tc.first];
        }
      }

      if (ttc.first < labels_size) {
        int tick_idx = ttc.first + fAdcDelay;
        if (tick_idx < labels_size) {
          labels_deposit[tick_idx] = max_deposit;
          labels_pdg[tick_idx] = max_pdg & type_pdg_mask;
        }
      }
    }

    for (auto const& drift_flags : wireToDriftToVtxFlags[widx]) {
      int drift = drift_flags.first, flags = drift_flags.second;
      if ((drift >= 0) && (drift < labels_size)) { labels_pdg[drift] |= flags; }
    }
    setWireEdepsAndLabels(labels_deposit, labels_pdg, widx);
  } // for each Wire

  return true;
}
// ------------------------------------------------------
*/
bool
nnet::TrainingDataAlg::findCrop(float max_e_cut,
                                unsigned int& w0,
                                unsigned int& w1,
                                unsigned int& d0,
                                unsigned int& d1) const
{
  if (debug_level>=0) std::cout << "PointIdAlg.cxx: findCrop " << std::endl;   
  if (fWireDriftEdep.empty() || fWireDriftEdep.front().empty()) return false;

  float max_cut = 0.25 * max_e_cut;

  w0 = 0;
  float cut = 0;
  while (w0 < fWireDriftEdep.size()) {
    for (auto const d : fWireDriftEdep[w0])
      cut += d;
    if (cut < max_cut)
      w0++;
    else
      break;
  }
  w1 = fWireDriftEdep.size() - 1;
  cut = 0;
  while (w1 > w0) {
    for (auto const d : fWireDriftEdep[w1])
      cut += d;
    if (cut < max_cut)
      w1--;
    else
      break;
  }
  w1++;

  d0 = 0;
  cut = 0;
  while (d0 < fWireDriftEdep.front().size()) {
    for (size_t i = w0; i < w1; ++i)
      cut += fWireDriftEdep[i][d0];
    if (cut < max_cut)
      d0++;
    else
      break;
  }
  d1 = fWireDriftEdep.front().size() - 1;
  cut = 0;
  while (d1 > d0) {
    for (size_t i = w0; i < w1; ++i)
      cut += fWireDriftEdep[i][d1];
    if (cut < max_cut)
      d1--;
    else
      break;
  }
  d1++;

  unsigned int margin = 20;
  if ((w1 - w0 > 8) && (d1 - d0 > 8)) {
    if (w0 < margin)
      w0 = 0;
    else
      w0 -= margin;

    if (w1 > fWireDriftEdep.size() - margin)
      w1 = fWireDriftEdep.size();
    else
      w1 += margin;

    if (d0 < margin)
      d0 = 0;
    else
      d0 -= margin;

    if (d1 > fWireDriftEdep.front().size() - margin)
      d1 = fWireDriftEdep.front().size();
    else
      d1 += margin;

    return true;
  }
  else
    return false;
}
// ------------------------------------------------------
