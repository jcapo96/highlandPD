#ifndef pionAnalysis_h
#define pionAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
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


class pionAnalysis: public baseAnalysis {
 public:
  pionAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~pionAnalysis(){}

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
  void InitializeBunch();
  void FillCategories();
  void DefineInputConverters();
  
  /// Returns the ToyBoxPD
  virtual const ToyBoxPD& box(Int_t isel=-1) const {return *static_cast<const ToyBoxPD*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}

  /// Method to add our own categories for color drawing
  void AddCustomCategories();

private:

public:

  // Needed to get the index of the counters from standardPDTree 
  Int_t seltrk_ndau;
  Int_t ntracks;  

  // Enum with unique indexes for output tree variables
  enum enumStandardMicroTrees_pionAnalysis{
    seltrk_truemom_tpc =     standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    seltrk_truebeta,
    seltrk_kinetic,
    
    seltrk_mom_muon_z0,
    seltrk_mom_prot_z0,
    seltrk_length_z0,
    seltrk_pos_z0,

    seltrk_broken,
    
    // seltrk_dedx_exp,

    seltrk_CNNscore,
    seltrk_chi2_prot,
    seltrk_chi2_ndf,
    seltrk_hit0_dqdx,
    seltrk_hit0_dedx,
    seltrk_truncLibo_dEdx,

    seltrk_dau_CNNscore,
    seltrk_dau_chi2_prot,
    seltrk_dau_chi2_ndf,
    seltrk_dau_vtxdistance,
    seltrk_dau_type,
    seltrk_dau_truncLibo_dEdx,
    seltrk_dau_hit0_dqdx,
    seltrk_dau_hit0_dedx,

    trk_pandora,
    enumStandardMicroTreesLast_pionAnalysis
  };

  enum enumConf_pionAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    Lifetime_syst,
    dQdxCalib_syst,
    Recombination_syst,
    enumConfLast_pionAnalysis
  };

  enum enumSyst_pionAnalysis{
    kdEdx=0,
    kLength,
    kLifetime,
    kdQdxCalib,
    kRecombination,
    kDerived,
    kBeam,    
    enumSystLast_pionAnalysis
  };
  
};

#endif
