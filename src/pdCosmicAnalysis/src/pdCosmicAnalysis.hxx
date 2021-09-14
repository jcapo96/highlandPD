#ifndef pdCosmicAnalysis_h
#define pdCosmicAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxKaon.hxx"
#include "standardPDTree.hxx"


/* This is an example of analysis in ProtoDUNE-SP detector 
   This example contains several elements:
     - A simple event selection of kaons decaying at rest from a beam of single kaons at 1GeV. Than means 
       that this example should be run over that particular sample
          /pnfs/dune/scratch/dunepro/v06_05_00/mergeana/gen_protoDune_kaon_1GeV_mono/13405798_0/anahist.root
     - A example of propagation of event variation systematic: dE/dx scaling
     - Filling of a root tree containing summary information about the selection and systematic mentioned above
     - A macro to produce few plots using the root file produced by the executable and the DrawingTools
 */


namespace pdCosmicAnalysisConstants{

  const UInt_t NMAXSAVEDCOSMICS        = 50;
  const UInt_t NMAXSAVEDTRUECANDIDATES = 10;
  const UInt_t NMAXSAVEDCANDIDATES     = 10;
}

class pdCosmicAnalysis: public baseAnalysis {
 public:
  pdCosmicAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~pdCosmicAnalysis(){}

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
  void InitializeBunch();
  void FillCategories();
  void DefineInputConverters();


  /// Get a casted AnaSpillC to AnaSpill from the InputManager
  AnaSpillPD& GetSpill(){return *static_cast<AnaSpillPD*>(&input().GetSpill());}
  
  /// Get a casted AnaBunchBB to AnaBunch from the InputManager
  AnaBunchPD& GetBunch(){return *static_cast<AnaBunchPD*>(&input().GetBunch());}

  /// Get a casted AnaEventC to AnaEvent 
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

  /// Method to add our own categories for color drawing
  void AddCustomCategories();

private:

public:

  enum enumConf_pdCosmicAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    Lifetime_syst,
    dQdxCalib_syst,
    Recombination_syst,
    dEdxCalib_syst,
    TrackEff_syst,
    enumConfLast_pdCosmicAnalysis
  };

  enum enumSyst_pdCosmicAnalysis{
    kLength=0,
    kLifetime,
    kdQdxCalib,
    kRecombination,
    kdEdxCalib,
    kBeam,
    kTrackEff,
    enumSystLast_pdCosmicAnalysis
  };  
};

#endif
