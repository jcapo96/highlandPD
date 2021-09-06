#ifndef pdExampleAnalysis_h
#define pdExampleAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"

/* This is an example of analysis in ProtoDUNE-SP detector 
   This example contains several elements:
     - A simple event selection of stopping protongs from a beam of single kaons at 1GeV. 
     - A example of propagation of event variation systematic: dE/dx scaling
     - Filling of a root tree containing summary information about the selection and systematic mentioned above
     - A macro to produce few plots using the root file produced by the executable and the DrawingTools
 */

namespace pdExampleAnalysisConstants{

  const UInt_t NMAXSAVEDPARTICLES      = 20;
  const UInt_t NMAXSAVEDDAUGHTERS      = 20;
  const UInt_t NMAXSAVEDGDAUGHTERS     = 50;
}


class pdExampleAnalysis: public baseAnalysis {
 public:
  pdExampleAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~pdExampleAnalysis(){}

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

  Int_t seltrk_ndau;
  Int_t ntracks;
  
public:

  enum enumStandardMicroTrees_pdExampleAnalysis{
    seltrk_csdarange_prot = standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    enumStandardMicroTreesLast_pdExampleAnalysis
  };

  enum enumConf_pdExampleAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    enumConfLast_pdExampleAnalysis
  };



};

#endif
