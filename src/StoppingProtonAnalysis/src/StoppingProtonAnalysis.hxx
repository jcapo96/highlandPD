#ifndef StoppingProtonAnalysis_h
#define StoppingProtonAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"

class StoppingProtonAnalysis: public baseAnalysis {
 public:
  StoppingProtonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~StoppingProtonAnalysis(){}

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

  using baseAnalysis::FillTruthTree;
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

  enum enumStandardMicroTrees_StoppingProtonAnalysis{
    seltrk_csdarange_prot = standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    seltrk_ndau,
    enumStandardMicroTreesLast_StoppingProtonAnalysis
  };

  enum enumConf_StoppingProtonAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    enumConfLast_StoppingProtonAnalysis
  };



};

#endif
