#ifndef pionMomentumAnalysis_h
#define pionMomentumAnalysis_h

#include "pdBaseAnalysis.hxx"
#include "standardPDTree.hxx"
#include "ToyBoxPD.hxx"

/// MicroTree + optional truth tree. When input is `mlskim` (skim ROOT), truth tree can record
/// e.g. the true PDG of the particle parenting secondaries (skim seed; expect 211 for pi+ skim).
/// For minitree input, truth-tree filling stays off unless you extend CheckFillTruthTreePD.
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

  enum enumPionMomentumTruthTreeExtra {
    pion_true_parent_of_secondaries_pdg = standardPDTree::enumStandardMicroTreesLast_standardPDTree + 1,
  };
};

#endif
