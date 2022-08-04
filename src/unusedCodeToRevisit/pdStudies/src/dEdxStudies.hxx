#ifndef dEdxStudies_h
#define dEdxStudies_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"

#include "pdCalorimetryUtils.hxx"
//#include "pdMVA.hxx"



/*

  This analysis algorithm is used to fill a tree with dEdx related variables computed in different ways, used for comparisons

 */


class dEdxStudies: public baseAnalysis {
 public:
  dEdxStudies(AnalysisAlgorithm* ana=NULL);
  virtual ~dEdxStudies(){}

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


    /// Create the appropriate event time from an Spill and a Bunch in that spill
  //  virtual AnaEventC* MakeEvent(){
  //    return new AnaEventPD(GetSpill(),GetBunch());
  //  }


  /// Get a casted AnaSpillC to AnaSpill from the InputManager
  AnaSpillPD& GetSpill(){return *static_cast<AnaSpillPD*>(&input().GetSpill());}


  void FilldEdxInfo(AnaParticlePD* part);
  
private:

  //  mvapid::MVAAlg fMVA;
  
public:

  // Needed to get the index of the counters from standardPDTree 
  Int_t seltrk_ndau;
  Int_t ntracks;  

  pdCalorimetryUtils fCalo;
  
  // Enum with unique indexes for output tree variables
  enum enumStandardMicroTrees_dEdxStudies{
    seltrk_nhits =     standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    seltrk_hit_dedx_0,
    seltrk_hit_dedx_1,
    seltrk_hit_dedx_2,
    seltrk_hit_dedx_3,
    seltrk_hit_dqdx_0,
    seltrk_hit_dqdx_3,
    seltrk_hit_dqdx_3_cal,
    seltrk_hit_dedx_0_cal,
    
    seltrk_hit_time,
    seltrk_hit_integral,
    seltrk_hit_amplitude,
    seltrk_hit_pitch3D,

    seltrk_dedx_fromdqdx,
    
    enumStandardMicroTreesLast_dEdxStudies
  };

  enum enumConf_dEdxStudies{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    enumConfLast_dEdxStudies
  };

  enum enumSyst_dEdxStudies{
    kdEdx=0,
    kLength,
    kBeam,
    enumSystLast_dEdxStudies
  };
  


};

#endif
