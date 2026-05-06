#ifndef pionMomentumAnalysis_h
#define pionMomentumAnalysis_h

#include "pdBaseAnalysis.hxx"
#include "pionMomentumEventDisplay.hxx"
#include "standardPDTree.hxx"
#include "ToyBoxPD.hxx"

/// MicroTree analysis for standard MiniTree inputs.
class pionMomentumAnalysis : public pdBaseAnalysis {
 public:
  explicit pionMomentumAnalysis(AnalysisAlgorithm* ana = NULL);
  virtual ~pionMomentumAnalysis();

  void DefineSelections();
  void DefineCorrections();
  void DefineSystematics();
  void DefineConfigurations();
  void DefineMicroTrees(bool addBase = true);
  void DefineTruthTree();

  void FillMicroTrees(bool addBase = true);
  void FillToyVarsInMicroTrees(bool addBase = true);

  bool CheckFillTruthTreePD(const AnaTrueParticlePD* part);

  using pdBaseAnalysis::FillTruthTree;
  void FillTruthTree(const AnaTrueParticlePD& part);

  bool Initialize();
  void Finalize();
  void FillCategories();
  void DefineInputConverters();

  virtual const ToyBoxPD& box(Int_t isel = -1) const {
    return *static_cast<const ToyBoxPD*>(&boxB(isel));
  }
  virtual AnaVertexB* GetVertex() const { return box().Vertex; }
  virtual AnaTrueVertexB* GetTrueVertex() const { return box().TrueVertex; }

 private:
  bool _ApplySCECorrection;
  bool _ApplySCESystematic;
  double _MCSRadiationLengthCm;
  double _MCSTargetSegmentLengthCm;
  double _MCSMinSegmentLengthCm;
  double _MCStheta0FloorRad;
  double _MCSMaxAbsDeltaThetaRad;
  int _MCSDropFirstNSegments;
  int _MCSDropLastNSegments;
  bool _MCSUseDetectorSigma;
  double _MCSDetectorSigmaFloorRad;
  double _MCSDetectorSigmaA;
  double _MCSDetectorSigmaC;
  double _MCSMomentumScanMaxGeV;
  int _TLEMinInteriorHits;
  int _TLESkipHitsFirst;
  int _TLESkipHitsLast;
  double _TLEDedxMinMeVcm;
  double _TLEDedxMaxMeVcm;
  double _TLEScanLmaxCm;
  double _TLEScanStepCm;
  double _TLEScanStepFineCm;
  double _TLELowPMomentumRefineGeV;
  double _TLEDedxPdfPathCm;
  double _BraggChi2PiMaxResidualRangeCm;
  double _BraggChi2PiSigmaMeVcm;
  int _BraggChi2PiMinHits;
  int _BraggChi2PiSkipHitsFirst;
  int _BraggChi2PiSkipHitsLast;
  double _BraggChi2PiDedxMinMeVcm;
  double _BraggChi2PiDedxMaxMeVcm;
  /// If true, run TLEFit and MCS momentum estimation only for true beam pions:
  /// abs(maintruepdg)==211 && maintrueid==beamtrueid.
  bool _MomentumRequireTrueBeamPion;
  bool _RunTLEFit;
  bool _RunMCS;
  bool _RunCSDA;
  bool _RunBeamDaughterTLETruncationScan;
  bool _RunBeamDaughterMCSTruncationScan;
  bool _EnableMomentumDiagnosticMultigraphs;
  bool _FreeRangeComputeDedxBiasDiagnostics;
  double _StoppingMaxTrueEndMomentumGeV;
  bool _CreateEventDisplay;
  pionMomentumEventDisplay* _EventDisplay;

 public:
  enum enumSyst_pionMomentumAnalysis {
    kSCEGeometric = 0,
    enumSystLast_pionMomentumAnalysis
  };

};

#endif
