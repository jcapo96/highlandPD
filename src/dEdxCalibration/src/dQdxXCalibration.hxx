#ifndef dQdxXCalibration_h
#define dQdxXCalibration_h

#include "TH3F.h"

#include "dQdxYZCalibration.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxdEdx.hxx"
#include "SpaceCharge.hxx"

class dQdxXCalibration: public baseAnalysis {
 public:
  dQdxXCalibration(AnalysisAlgorithm* ana=NULL);
  virtual ~dQdxXCalibration(){}

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
  void FillHistograms();
  void FillToyHistograms();

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

  int _SelectedTracks;
  int _MaxTracks;

  bool _SaveAna;
  bool _SaveToy;

  bool _ApplySCEPositionCorrection;
  bool _ApplySCEPitchCorrection;
  bool _ApplyLifetimeCorrection;

  bool _ApplySCESystematic;
  bool _ApplyLifetimeSystematic;

  int _MichelRemovingTree;

  static const int nbinsx = (XMAX-XMIN)/STEP;
  TH1F* h_global_x;
  TH1F* h_local_x[nbinsx];
  TH2F* h_global_x_toy;
  TH2F* h_local_x_toy[nbinsx];

  TH2F* yz_correction;
  TH3F* yz_correction_toy;

  SpaceCharge* _sce;

  SpaceCharge* GetSCE();
  
public:

  enum enumSyst_dQdxXCalibration{
    kSCE_syst=0,
    kLifetime_syst,
    enumSystLast_dQdxXCalibration
  };

  enum enumCorr_dQdxXCalibration{
    kSCE_corr=0,
    kSCEPosition_corr,
    kSCEPitch_corr,
    kLifetime_corr,
    enumCorrLast_dQdxXCalibration
  };

  enum enumStandardMicroTrees_dQdxXCalibration{
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
    enumStandardMicroTreesLast_dQdxXCalibration
  };

  

  
};

#endif
