#ifndef neutralKaonAnalysis_h
#define neutralKaonAnalysis_h

#include "pdBaseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"
#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysisUtils.hxx"

namespace neutralKaonAnalysisConstants{

  const UInt_t NMAXSAVEDCANDIDATES     = 10;
}

class neutralKaonAnalysis: public pdBaseAnalysis {
 public:
 neutralKaonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~neutralKaonAnalysis(){
  }

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
  void FillTruthTree(const AnaTrueParticlePD& part);
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
  // SCE correction and systematic parameters
  bool _ApplySCECorrection;
  bool _ApplySCESystematic;


public:

  enum enumStandardMicroTrees_neutralKaonAnalysis{
    seltrk_csdarange_prot = neutralKaonTree::enumNeutralKaonMicroTreesLast,
    seltrk_ndau,
    seltrk_dau_ndau,
    seltrk_truthdau_ndau,
    seltrk_dau_trueparentpdg,

    // Debugging variables
    nAllParticles,

    enumStandardMicroTreesLast_neutralKaonAnalysis
  };

  enum enumConf_neutralKaonAnalysis{
    detmass_syst=pdBaseAnalysis::enumConfLast_baseAnalysis+1,
    dedx_syst,
    enumConfLast_neutralKaonAnalysis
  };

  enum enumSyst_neutralKaonAnalysis{
    kSCEGeometric=0,
    enumSystLast_neutralKaonAnalysis
  };


  // Check if event contains signal vertices (K0 -> pi+ pi-)
  bool EventContainsSignalVertices(const AnaEventB& event) const;

};

#endif
