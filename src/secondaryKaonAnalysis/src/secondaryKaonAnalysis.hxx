#ifndef secondaryKaonAnalysis_h
#define secondaryKaonAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxKaon.hxx"
#include "standardPDTree.hxx"
#include "kaonTree.hxx"
#include "kaonAnalysisUtils.hxx"
#include "kaonDataClasses.hxx"


namespace secondaryKaonAnalysisConstants{

  const UInt_t NMAXSAVEDPARTICLES      = 20;
  const UInt_t NMAXSAVEDDAUGHTERS      = 20;
  const UInt_t NMAXSAVEDGDAUGHTERS     = 5;
  const UInt_t NMAXSAVEDGGDAUGHTERS    = 5;
  const UInt_t NMAXSAVEDHITSDAUGHTERS  = 6000;
  const UInt_t NMAXSAVEDTRUECANDIDATES = 10;
  const UInt_t NMAXSAVEDCANDIDATES     = 10;
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

  bool CheckFillTruthTree(const AnaTrueVertex& vtx);

  using baseAnalysis::FillTruthTree;
  void FillTruthTree(const AnaTrueVertex& vtx);
  //--------------------

  bool Initialize();
  bool InitializeSpill(); //temporary solution for filling truthtree
  void FinalizeSpill(); //temporary solution for filling truthtree
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
  
  bool _AllSelection;
  bool _BeamSelection;
  bool _CosmicSelection;

public:

  enum enumConf_secondaryKaonAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    Lifetime_syst,
    dQdxCalib_syst,
    Recombination_syst,
    dEdxCalib_syst,
    TrackEff_syst,
    enumConfLast_secondaryKaonAnalysis
  };

  enum enumSyst_secondaryKaonAnalysis{
    kLength=0,
    kLifetime,
    kdQdxCalib,
    kRecombination,
    kdEdxCalib,
    kBeam,
    kTrackEff,
    enumSystLast_secondaryKaonAnalysis
  };
  
};

#endif
