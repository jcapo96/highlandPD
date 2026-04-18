#ifndef pionMomentumSkimAnalysis_h
#define pionMomentumSkimAnalysis_h

#include "pdBaseAnalysis.hxx"
#include "pionMomentumMlSkim.hxx"
#include "ToyBoxPD.hxx"

/// Skim-only driver: writes `mlskim` via pionMomentumMlSkim. Run with AnalysisLoop `-r` so no main microTree
/// file is produced (`-o` is still required by the CLI but not opened when `-r` is set).
class pionMomentumSkimAnalysis : public pdBaseAnalysis {
 public:
  explicit pionMomentumSkimAnalysis(AnalysisAlgorithm* ana = NULL);
  ~pionMomentumSkimAnalysis();

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
  void FinalizeBunch() override;
  void FillCategories();
  void DefineInputConverters();

  virtual const ToyBoxPD& box(Int_t isel = -1) const {
    return *static_cast<const ToyBoxPD*>(&boxB(isel));
  }
  virtual AnaVertexB* GetVertex() const { return box().Vertex; }
  virtual AnaTrueVertexB* GetTrueVertex() const { return box().TrueVertex; }

 private:
  pionMomentumMlSkim _skim;
};

#endif
