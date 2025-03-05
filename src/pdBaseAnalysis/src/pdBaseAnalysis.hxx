#ifndef pdBaseAnalysis_h
#define pdBaseAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"

namespace pdBaseAnalysisConstants{

  //const UInt_t NMAXSAVEDCANDIDATES     = 20;
}

class pdBaseAnalysis: public baseAnalysis {
 public:
  pdBaseAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~pdBaseAnalysis(){}

  //---- These are mandatory functions. Dummy implementation for the moment
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

  bool FinalizeConfiguration();
  void FillTruthTree();
  void FillTruthTree(const AnaTrueVertex& vtx){return;}
  virtual void FillTruthTree(const AnaTrueParticlePD& part);
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
  
  /// Returns the ToyBoxPD
  virtual const ToyBoxPD& box(Int_t isel=-1) const {return *static_cast<const ToyBoxPD*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaParticlePD* GetParticle() const{return box().MainTrack;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueParticlePD* GetTrueParticle() const {return box().TrueMainTrack;}


private:

public:

  enum enumConf_pdBaseAnalysis{
    enumConfFirst_pdBaseAnalysis=baseAnalysis::enumConfLast_baseAnalysis+1,    
    enumConfLast_pdBaseAnalysis
  };

  enum enumSyst_pdBaseAnalysis{
    enumSystFirst_pdBaseAnalysis=0,
    enumSystLast_pdBaseAnalysis
  };
   
};

#endif
