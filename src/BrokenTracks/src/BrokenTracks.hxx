#ifndef BrokenTracks_h
#define BrokenTracks_h

#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxdEdx.hxx"

class BrokenTracks: public baseAnalysis {
 public:
  BrokenTracks(AnalysisAlgorithm* ana=NULL);
  virtual ~BrokenTracks(){}

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

  using baseAnalysis::FillTruthTree;
  void FillTruthTree(const AnaTrueVertex& vtx);
  //--------------------

  bool Initialize();
  void FillCategories();
  void DefineInputConverters();

  void AddBrokenCategory();
  void FillBrokenCategory(AnaParticlePD* part);

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
  
  /// Returns the ToyBox
  virtual const ToyBoxdEdx& box(Int_t isel=-1) const {return *static_cast<const ToyBoxdEdx*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}
  
private:

public:

  enum enumSyst_BrokenTracks{
    kSCE_syst=0,
    kLifetime_syst,
    enumSystLast_BrokenTracks
  };

  enum enumCorr_BrokenTracks{
    kSCE_corr=0,
    kSCEPosition_corr,
    kSCEPitch_corr,
    kLifetime_corr,
    enumCorrLast_BrokenTracks
  };

};

#endif
