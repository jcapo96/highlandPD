#ifndef dQdxYZCalibration_h
#define dQdxYZCalibration_h

#include "TH2F.h"

#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxdEdx.hxx"

const int NTOYS = 100;
const int ZMIN  = 0;
const int ZMAX  = 695;
const int YMIN  = 0;
const int YMAX  = 600;
const int XMIN  = -360;
const int XMAX  = 0;
const int STEP  = 5; 

class dQdxYZCalibration: public baseAnalysis {
 public:
  dQdxYZCalibration(AnalysisAlgorithm* ana=NULL);
  virtual ~dQdxYZCalibration(){}

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

  void Finalize();

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

  bool _SaveAna;
  bool _SaveToy;

  bool _ApplySCEPositionCorrection;
  bool _ApplySCEPitchCorrection;
  bool _ApplyLifetimeCorrection;

  bool _ApplySCESystematic;
  bool _ApplyLifetimeSystematic;

  int _MichelRemovingTree;

  static const int nbinsy = (YMAX-YMIN)/STEP;
  static const int nbinsz = (ZMAX-ZMIN)/STEP;
  TH2F* h_global_yz;
  TH2F* h_local_yz[nbinsz][nbinsy];
  
public:

  enum enumSyst_dQdxYZCalibration{
    kSCE_syst=0,
    kLifetime_syst,
    enumSystLast_dQdxYZCalibration
  };

  enum enumCorr_dQdxYZCalibration{
    kSCE_corr=0,
    kSCEPosition_corr,
    kSCEPitch_corr,
    kLifetime_corr,
    enumCorrLast_dQdxYZCalibration
  };

  enum enumStandardMicroTrees_dQdxYZCalibration{
    ntracks = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    track_hit_x,
    track_hit_y,
    track_hit_z,
    track_hit_dqdx,
    toy_track_hit_x,
    toy_track_hit_y,
    toy_track_hit_z,
    toy_track_hit_dqdx,
    hit_x,
    hit_y,
    hit_z,
    hit_dqdx,
    toy_hit_x,
    toy_hit_y,
    toy_hit_z,
    toy_hit_dqdx,
    enumStandardMicroTreesLast_dQdxYZCalibration
  };

  

  
};

#endif
