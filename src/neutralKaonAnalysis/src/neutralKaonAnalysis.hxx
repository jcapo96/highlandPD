#ifndef neutralKaonAnalysis_h
#define neutralKaonAnalysis_h

#include "pdBaseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"

class neutralKaonAnalysis: public pdBaseAnalysis {
 public:
 neutralKaonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~neutralKaonAnalysis(){}

  //---- These are mandatory functions
  void DefineSelections();
  void DefineCorrections();
  void DefineSystematics();
  void DefineConfigurations();
  void DefineMicroTrees(bool addBase=true);
  void DefineTruthTree();

  void FillMicroTrees(bool addBase=true);
  void FillToyVarsInMicroTrees(bool addBase=true);

  bool CheckFillTruthTree(const AnaTrueVertex& vtx);
  bool CheckFillTruthTreePD(const AnaTrueParticlePD* part);

  using pdBaseAnalysis::FillTruthTree;
  void FillTruthTree(const AnaTrueVertex& vtx);
  //--------------------

  bool Initialize();
  void FillCategories();
  void DefineInputConverters();

  /// Returns the ToyBoxPD
  virtual const ToyBoxPD& box(Int_t isel=-1) const {return *static_cast<const ToyBoxPD*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}


private:

public:

  enum enumStandardMicroTrees_neutralKaonAnalysis{
    seltrk_csdarange_prot = standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    seltrk_ndau,
    enumStandardMicroTreesLast_neutralKaonAnalysis
  };

  enum enumConf_neutralKaonAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,
    dedx_syst,
    enumConfLast_neutralKaonAnalysis
  };



};

#endif
