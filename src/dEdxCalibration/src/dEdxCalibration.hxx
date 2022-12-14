#ifndef dEdxCalibration_h
#define dEdxCalibration_h

#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxdEdx.hxx"

class dEdxCalibration: public baseAnalysis {
 public:
  dEdxCalibration(AnalysisAlgorithm* ana=NULL);
  virtual ~dEdxCalibration(){}

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

  enum enumStandardMicroTrees_dEdxCalibration{
    ntracks = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    track_hit_x,
    track_hit_y,
    track_hit_z,
    track_hit_dqdx,
    // toy_track_hit_x,
    // toy_track_hit_y,
    // toy_track_hit_z,
    // toy_track_hit_dqdx,
    enumStandardMicroTreesLast_dEdxCalibration
  };

  
};

#endif
