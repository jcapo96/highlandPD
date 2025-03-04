#ifndef secondaryKaonAnalysis_h
#define secondaryKaonAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxKaon.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "kaonTree.hxx"
#include "kaonAnalysisUtils.hxx"

#include "TH1F.h"

namespace secondaryKaonAnalysisConstants{

  const UInt_t NMAXSAVEDCANDIDATES     = 20;
}

class secondaryKaonAnalysis: public baseAnalysis {
 public:
  secondaryKaonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~secondaryKaonAnalysis(){}

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

  bool _xs_selection;
  std::string _selection_name;
  
  bool _UseDetailedSelection;

  bool _FillBeamInstrumentationInfo;
  bool _FillBeamParticleInfo;
  bool _FillBeamParticleHitsInfo;
  bool _FillBeamParticleDaughtersInfo;
  bool _FillKaonCandidatesInfo;
  bool _FillBestKaonCandidateInfo;
  bool _FillToyVariablesInfo;

  bool _ApplydQdxSystematic;
  bool _ApplyRecombinationSystematic;
  bool _ApplySCESystematic;
  bool _ApplyBrokenTracksSystematic;
  bool _ApplyBeamPIDEfficiencySystematic;
  bool _ApplyBeamPartWeightSystematic;
  bool _ApplyBeamMomWeightSystematic;
  bool _ApplyProtonBackgroundWeightSystematic;

  bool _ApplySCECorrection;
  bool _ApplydEdxCorrection;
  bool _ApplyRecombinationCorrection;

public:

  enum enumConf_secondaryKaonAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    enumConfLast_secondaryKaonAnalysis
  };

  enum enumSyst_secondaryKaonAnalysis{
    kdQdxCalibration=0,
    kRecombination,
    kSCEGeometric,
    kBrokenTracks,
    kBeamPIDEfficiency,
    kBeamPartWeight,
    kBeamMomWeight,
    kProtonBackgroundWeight,
    enumSystLast_secondaryKaonAnalysis
  };
  
  enum enumKaonSystVarMicroTrees{

    //true kaon candidates info
    bestcandidate_hit_resrange_toy = kaonTree::enumKaonMicroTreesLast+1,
    bestcandidate_hit_dedx_toy,
    bestcandidate_dau_hit_resrange_toy,
    bestcandidate_dau_hit_dedx_toy,
    bestcandidate_chi2_kaon_perndf_toy,
    bestcandidate_dau_chi2_muon_perndf_toy,
    bestcandidate_Zstart_toy,
    bestcandidate_Zend_toy,
    bestcandidate_calE_toy,
    bestcandidate_length_toy,
    bestcandidate_dau_calE_toy,
  };
  
};

#endif
