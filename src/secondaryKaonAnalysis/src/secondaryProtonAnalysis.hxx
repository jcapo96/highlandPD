#ifndef secondaryProtonAnalysis_h
#define secondaryProtonAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxKaon.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "kaonTree.hxx"
#include "kaonAnalysisUtils.hxx"

#include "TH1F.h"

class secondaryProtonAnalysis: public baseAnalysis {
 public:
  secondaryProtonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~secondaryProtonAnalysis(){}

  //---- These are mandatory functions
  void DefineSelections();
  void DefineCorrections();
  void DefineSystematics();
  void DefineConfigurations();
  void DefineMicroTrees(bool addBase=true);
  void DefineTruthTree();

  void FillMicroTrees(bool addBase=true);
  void FillToyVarsInMicroTrees(bool addBase=true);

  bool CheckFillTruthTree(const AnaTrueVertex& vtx){return true;}
  bool CheckFillTruthTree(const AnaTrueParticlePD* part);

  //using baseAnalysis::FillTruthTree;
  bool FinalizeConfiguration();
  void FillTruthTree();
  void FillTruthTree(const AnaTrueVertex& vtx){return;}
  void FillTruthTree(const AnaTrueParticlePD& part);
  //--------------------

  bool Initialize();
  void FillCategories();
  void DefineInputConverters();


  /// Get a casted AnaSpillC to AnaSpillPD from the InputManager
  AnaSpillPD& GetSpill(){return *static_cast<AnaSpillPD*>(&input().GetSpill());}
  
  /// Get a casted AnaBunchBB to AnaBunchPD from the InputManager
  AnaBunchPD& GetBunch(){return *static_cast<AnaBunchPD*>(&input().GetBunch());}

  /// Get a casted AnaEventC to AnaEventPD 
  AnaEventPD& GetEvent(){return *static_cast<AnaEventPD*>(_event);}
  
  /// Create the appropriate event time from an Spill and a Bunch in that spill
  virtual AnaEventC* MakeEvent(){
    return new AnaEventPD(GetSpill(),GetBunch());
  }
  
  /// Returns the ToyBoxKaon
  virtual const ToyBoxKaon& box(Int_t isel=-1) const {return *static_cast<const ToyBoxKaon*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxKaon
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxKaon
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}


private:

  enum enumSecondaryProtonVarMicroTrees{
    candidates_allhits = kaonTree::enumKaonMicroTreesLast+1,
    candidates_hitvector_dedx,
    candidates_hitvector_resrange,
    candidates_hitvector_thetaYZ,
    candidates_hitvector_thetaXZ,
    candidates_hitvector_tracklength,
  };

public:

};

#endif
