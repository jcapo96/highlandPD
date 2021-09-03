/*
   converter for trees created by the official PDSPAnalyzer module

   M. García June 2021
*/

#ifndef PDSPAnalyzerTreeConverter_h
#define PDSPAnalyzerTreeConverter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <set>
#include "InputConverter.hxx"
//#include "GeometryManager.hxx"

#include "pdDataClasses.hxx"

using namespace std;

class PDSPAnalyzerTreeConverter: public InputConverter{

 public:

  PDSPAnalyzerTreeConverter();
  virtual ~PDSPAnalyzerTreeConverter();

  virtual bool Initialize();
  virtual Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill);
  Int_t GetEvent(Long64_t& entry, AnaEventC*& event){(void)entry;(void)event; return 0;}

  /// Record the POT for the current spill, based on information in the AnaBeam
  /// member of the current AnaSpill.
  void IncrementPOTBySpill(){return;}

  virtual Int_t ReadEntries(Long64_t& entry);
  virtual bool AddFileToTChain(const std::string& inputString);
  void InitializeVariables();
  void SetBranchAddresses();

  //----------------
  virtual AnaSpillB*          MakeSpill()       { return new AnaSpill(); }
  virtual AnaBunch*           MakeBunch()       { return new AnaBunch(); }
  virtual AnaBeamPD*          MakeBeam()        { return new AnaBeamPD(); }
  virtual AnaDataQualityB*    MakeDataQuality() { return new AnaDataQuality(); }
  virtual AnaEventInfoB*      MakeEventInfo()   { return new AnaEventInfo(); }
  virtual AnaTrigger*         MakeTrigger()     { return new AnaTrigger(); }

  virtual AnaTrueParticlePD*  MakeTrueParticle(){ return new AnaTrueParticlePD(); }
  virtual AnaTrueVertex*      MakeTrueVertex()  { return new AnaTrueVertex(); }
  virtual AnaParticlePD*      MakeParticle()    { return new AnaParticlePD(); }

  // ----------------------------

  virtual void FillInfo(AnaSpill* spill);

  virtual void FillDQInfo(AnaDataQuality* dq);

  virtual void FillTrueInfo(AnaSpill* spill);
  virtual void FillTrueBeamTrueParticleInfo(AnaTrueParticlePD* truePart);
  virtual void FillTrueBeamDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart, AnaTrueParticlePD* parentPart);
  virtual void FillTrueBeamGrandDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart);

  virtual void FillBeamInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBeamPD* beam);

  virtual void FillBunchInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBunch* bunch, AnaBeamPD* beam);

  virtual void FillBeamParticleInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaParticlePD* part, AnaBeamPD* beam);
  virtual void FillDaughterPFPInfo(Int_t itrk, AnaParticlePD* part);
  virtual void FillDaughterParticleTrackInfo(std::vector<AnaTrueParticleB*>& trueParticles, Int_t itrk, AnaParticlePD* part);
  virtual void FillDaughterParticleShowerInfo(std::vector<AnaTrueParticleB*>& trueParticles, Int_t itrk, AnaParticlePD* part);

  virtual void FillBeamTrueParticleInfo(AnaTrueParticlePD* truePart);
  virtual void FillDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart);
  
  AnaTrueObjectC* FindTrueParticle(Int_t itrk, std::vector<AnaTrueParticleB*>& trueParticles);

protected:

  AnaSpill* _spill;
  
  std::string _previousFile;
  Int_t _previousRunID;
  Int_t _previousSubrunID;
  Int_t _previousRefEventID;

  bool _byHits;
  
 protected:

  // TChains   
  TChain *eventsTree;
  TChain *FileIndexTree;

  Int_t Entries; 
  Int_t Counter; 

  Bool_t _isMC;
  std::string _softwareVersion;

  // Declaration of leaf types
  // Meta data
  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Int_t           MC;
  // True beam info
    Int_t           true_beam_PDG;
  Double_t        true_beam_mass;
  Int_t           true_beam_ID;
  string          *true_beam_endProcess;
  Double_t        true_beam_startX;
  Double_t        true_beam_startY;
  Double_t        true_beam_startZ;
  Double_t        true_beam_endX;
  Double_t        true_beam_endY;
  Double_t        true_beam_endZ;
  Double_t        true_beam_startP;
  Double_t        true_beam_startPx;
  Double_t        true_beam_startPy;
  Double_t        true_beam_startPz;
  Double_t        true_beam_endP;
  Double_t        true_beam_endPx;
  Double_t        true_beam_endPy;
  Double_t        true_beam_endPz;
  Double_t        true_beam_last_len;
  Double_t        true_beam_startDirX;
  Double_t        true_beam_startDirY;
  Double_t        true_beam_startDirZ;
  Int_t           true_beam_nElasticScatters;
  vector<double>  *true_beam_elastic_costheta;
  vector<double>  *true_beam_elastic_X;
  vector<double>  *true_beam_elastic_Y;
  vector<double>  *true_beam_elastic_Z;
  vector<double>  *true_beam_elastic_deltaE;
  vector<double>  *true_beam_elastic_IDE_edep;
  Double_t        true_beam_IDE_totalDep;
  Int_t           true_daughter_nPi0;
  Int_t           true_daughter_nPiPlus;
  Int_t           true_daughter_nProton;
  Int_t           true_daughter_nNeutron;
  Int_t           true_daughter_nPiMinus;
  Int_t           true_daughter_nNucleus;
  vector<string>  *true_beam_processes;
  vector<double>  *true_beam_incidentEnergies;
  Double_t        true_beam_interactingEnergy;
  vector<int>     *true_beam_slices;
  vector<double>  *true_beam_slices_deltaE;
  vector<double>  *true_beam_traj_X;
  vector<double>  *true_beam_traj_Y;
  vector<double>  *true_beam_traj_Z;
  vector<double>  *true_beam_traj_KE;
  // Reco beam info
  Int_t           reco_beam_type;
  Int_t           reco_beam_trackID;
  Double_t        reco_beam_startX;
  Double_t        reco_beam_startY;
  Double_t        reco_beam_startZ;
  Double_t        reco_beam_endX;
  Double_t        reco_beam_endY;
  Double_t        reco_beam_endZ;
  Double_t        reco_beam_trackDirX;
  Double_t        reco_beam_trackDirY;
  Double_t        reco_beam_trackDirZ;
  Double_t        reco_beam_trackEndDirX;
  Double_t        reco_beam_trackEndDirY;
  Double_t        reco_beam_trackEndDirZ;
  Bool_t          reco_beam_flipped;
  Double_t        reco_beam_len;
  Double_t        reco_beam_alt_len;
  Double_t        reco_beam_momByRange_proton;
  Double_t        reco_beam_momByRange_muon;
  Double_t        reco_beam_momByRange_alt_proton;
  Double_t        reco_beam_momByRange_alt_muon;
  Double_t        reco_beam_calo_startX;
  Double_t        reco_beam_calo_startY;
  Double_t        reco_beam_calo_startZ;
  Double_t        reco_beam_calo_endX;
  Double_t        reco_beam_calo_endY;
  Double_t        reco_beam_calo_endZ;
  Double_t        reco_beam_calo_startDirX;
  Double_t        reco_beam_calo_startDirY;
  Double_t        reco_beam_calo_startDirZ;
  Double_t        reco_beam_calo_endDirX;
  Double_t        reco_beam_calo_endDirY;
  Double_t        reco_beam_calo_endDirZ;
  Double_t        reco_beam_vertex_michel_score;
  Int_t           reco_beam_vertex_nHits;
  vector<double>  *reco_beam_dQdX_SCE;
  vector<double>  *reco_beam_dQ;
  vector<double>  *reco_beam_dEdX_SCE;
  vector<double>  *reco_beam_calibrated_dEdX_SCE;
  vector<double>  *reco_beam_resRange_SCE;
  vector<double>  *reco_beam_TrkPitch_SCE;
  vector<double>  *reco_beam_dQdX_NoSCE;
  vector<double>  *reco_beam_dEdX_NoSCE;
  vector<double>  *reco_beam_calibrated_dEdX_NoSCE;
  vector<double>  *reco_beam_resRange_NoSCE;
  vector<double>  *reco_beam_TrkPitch_NoSCE;
  vector<double>  *reco_beam_calo_wire;
  vector<double>  *reco_beam_calo_wire_z;
  vector<double>  *reco_beam_calo_tick;
  vector<int>     *reco_beam_calo_TPC;
  Bool_t          reco_beam_passes_beam_cuts;
  Double_t        reco_beam_Chi2_proton;
  Int_t           reco_beam_Chi2_ndof;
  Int_t           reco_beam_PFP_ID;
  Int_t           reco_beam_PFP_nHits;
  Double_t        reco_beam_PFP_trackScore;
  Double_t        reco_beam_PFP_emScore;
  Double_t        reco_beam_PFP_michelScore;
  Double_t        reco_beam_PFP_trackScore_collection;
  Double_t        reco_beam_PFP_emScore_collection;
  Double_t        reco_beam_PFP_michelScore_collection;
  Int_t           reco_beam_allTrack_ID;
  Bool_t          reco_beam_allTrack_beam_cuts;
  Bool_t          reco_beam_allTrack_flipped;
  Double_t        reco_beam_allTrack_len;
  Double_t        reco_beam_allTrack_startX;
  Double_t        reco_beam_allTrack_startY;
  Double_t        reco_beam_allTrack_startZ;
  Double_t        reco_beam_allTrack_endX;
  Double_t        reco_beam_allTrack_endY;
  Double_t        reco_beam_allTrack_endZ;
  Double_t        reco_beam_allTrack_trackDirX;
  Double_t        reco_beam_allTrack_trackDirY;
  Double_t        reco_beam_allTrack_trackDirZ;
  Double_t        reco_beam_allTrack_trackEndDirX;
  Double_t        reco_beam_allTrack_trackEndDirY;
  Double_t        reco_beam_allTrack_trackEndDirZ;
  vector<double>  *reco_beam_allTrack_resRange;
  vector<double>  *reco_beam_allTrack_calibrated_dEdX;
  Double_t        reco_beam_allTrack_Chi2_proton;
  Int_t           reco_beam_allTrack_Chi2_ndof;
  string          *reco_beam_true_byHits_endProcess;
  string          *reco_beam_true_byHits_process;
  Int_t           reco_beam_true_byHits_origin;
  Int_t           reco_beam_true_byHits_PDG;
  Int_t           reco_beam_true_byHits_ID;
  Bool_t          reco_beam_true_byHits_matched;
  Double_t        reco_beam_true_byHits_purity;
  Double_t        reco_beam_true_byHits_startE;
  Double_t        reco_beam_true_byHits_startP;
  Double_t        reco_beam_true_byHits_startPx;
  Double_t        reco_beam_true_byHits_startPy;
  Double_t        reco_beam_true_byHits_startPz;
  Double_t        reco_beam_true_byHits_endE;
  Double_t        reco_beam_true_byHits_endP;
  Double_t        reco_beam_true_byHits_endPx;
  Double_t        reco_beam_true_byHits_endPy;
  Double_t        reco_beam_true_byHits_endPz;
  Double_t        reco_beam_true_byE_startE;
  Double_t        reco_beam_true_byE_startP;
  Double_t        reco_beam_true_byE_startPx;
  Double_t        reco_beam_true_byE_startPy;
  Double_t        reco_beam_true_byE_startPz;
  Double_t        reco_beam_true_byE_endE;
  Double_t        reco_beam_true_byE_endP;
  Double_t        reco_beam_true_byE_endPx;
  Double_t        reco_beam_true_byE_endPy;
  Double_t        reco_beam_true_byE_endPz;
  string          *reco_beam_true_byE_endProcess;
  string          *reco_beam_true_byE_process;
  Int_t           reco_beam_true_byE_origin;
  Int_t           reco_beam_true_byE_PDG;
  Int_t           reco_beam_true_byE_ID;
  Bool_t          reco_beam_true_byE_matched;
  vector<double>  *reco_beam_incidentEnergies;
  Double_t        reco_beam_interactingEnergy;
  // Reco info for all tracks
  vector<double>          *reco_track_startX;
  vector<double>          *reco_track_startY;
  vector<double>          *reco_track_startZ;
  vector<double>          *reco_track_endX;
  vector<double>          *reco_track_endY;
  vector<double>          *reco_track_endZ;
  vector<double>          *reco_track_michel_score;
  vector<int>             *reco_track_nHits;
  vector<int>             *reco_track_ID;
  // True daughter info
  vector<int>     *true_beam_daughter_ID;
  vector<int>     *true_beam_daughter_PDG;
  vector<double>  *true_beam_daughter_len;
  vector<double>  *true_beam_daughter_startX;
  vector<double>  *true_beam_daughter_startY;
  vector<double>  *true_beam_daughter_startZ;
  vector<double>  *true_beam_daughter_endX;
  vector<double>  *true_beam_daughter_endY;
  vector<double>  *true_beam_daughter_endZ;
  vector<double>  *true_beam_daughter_startP;
  vector<double>  *true_beam_daughter_startPx;
  vector<double>  *true_beam_daughter_startPy;
  vector<double>  *true_beam_daughter_startPz;
  vector<string>  *true_beam_daughter_Process;
  vector<string>  *true_beam_daughter_endProcess;
  vector<int>     *true_beam_grand_daughter_ID;
  vector<int>     *true_beam_grand_daughter_parID;
  vector<int>     *true_beam_grand_daughter_PDG;
  vector<string>  *true_beam_grand_daughter_Process;
  vector<string>  *true_beam_grand_daughter_endProcess;
  vector<int>     *true_beam_Pi0_decay_ID;
  vector<int>     *true_beam_Pi0_decay_parID;
  vector<int>     *true_beam_Pi0_decay_PDG;
  vector<double>  *true_beam_Pi0_decay_startX;
  vector<double>  *true_beam_Pi0_decay_startY;
  vector<double>  *true_beam_Pi0_decay_startZ;
  vector<double>  *true_beam_Pi0_decay_startP;
  vector<double>  *true_beam_Pi0_decay_startPx;
  vector<double>  *true_beam_Pi0_decay_startPy;
  vector<double>  *true_beam_Pi0_decay_startPz;
  vector<double>  *true_beam_Pi0_decay_len;
  //reco daughter info
  vector<int>             *reco_daughter_PFP_ID;
  vector<int>             *reco_daughter_PFP_nHits;
  vector<double>          *reco_daughter_PFP_trackScore;
  vector<double>          *reco_daughter_PFP_emScore;
  vector<double>          *reco_daughter_PFP_michelScore;
  vector<int>             *reco_daughter_PFP_nHits_collection;
  vector<double>          *reco_daughter_PFP_trackScore_collection;
  vector<double>          *reco_daughter_PFP_emScore_collection;
  vector<double>          *reco_daughter_PFP_michelScore_collection;
  vector<int>             *reco_daughter_allTrack_ID;
  vector<double>          *reco_daughter_allTrack_Theta;
  vector<double>          *reco_daughter_allTrack_Phi;
  vector<double>          *reco_daughter_allTrack_len;
  vector<double>          *reco_daughter_allTrack_alt_len;
  vector<double>          *reco_daughter_allTrack_momByRange_proton;
  vector<double>          *reco_daughter_allTrack_momByRange_muon;
  vector<double>          *reco_daughter_allTrack_momByRange_alt_proton;
  vector<double>          *reco_daughter_allTrack_momByRange_alt_muon;
  vector<double>          *reco_daughter_allTrack_startX;
  vector<double>          *reco_daughter_allTrack_startY;
  vector<double>          *reco_daughter_allTrack_startZ;
  vector<double>          *reco_daughter_allTrack_endX;
  vector<double>          *reco_daughter_allTrack_endY;
  vector<double>          *reco_daughter_allTrack_endZ;
  vector<double>          *reco_daughter_allTrack_vertex_michel_score;
  vector<int>             *reco_daughter_allTrack_vertex_nHits;
  vector<vector<double> > *reco_daughter_allTrack_dQdX_SCE;
  vector<vector<double> > *reco_daughter_allTrack_dEdX_SCE;
  vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE;
  vector<vector<double> > *reco_daughter_allTrack_resRange_SCE;
  vector<double>          *reco_daughter_allTrack_Chi2_proton;
  vector<int>             *reco_daughter_allTrack_Chi2_ndof;
  vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE_plane0;
  vector<vector<double> > *reco_daughter_allTrack_resRange_plane0;
  vector<double>          *reco_daughter_allTrack_Chi2_proton_plane0;
  vector<int>             *reco_daughter_allTrack_Chi2_ndof_plane0;
  vector<vector<double> > *reco_daughter_allTrack_calibrated_dEdX_SCE_plane1;
  vector<vector<double> > *reco_daughter_allTrack_resRange_plane1;
  vector<double>          *reco_daughter_allTrack_Chi2_proton_plane1;
  vector<int>             *reco_daughter_allTrack_Chi2_ndof_plane1;
  vector<int>             *reco_daughter_allShower_ID;
  vector<double>          *reco_daughter_allShower_len;
  vector<double>          *reco_daughter_allShower_startX;
  vector<double>          *reco_daughter_allShower_startY;
  vector<double>          *reco_daughter_allShower_startZ;
  vector<double>          *reco_daughter_allShower_dirX;
  vector<double>          *reco_daughter_allShower_dirY;
  vector<double>          *reco_daughter_allShower_dirZ;
  vector<double>          *reco_daughter_allShower_energy;
  vector<int>             *reco_daughter_PFP_true_byHits_PDG;
  vector<int>             *reco_daughter_PFP_true_byHits_ID;
  vector<int>             *reco_daughter_PFP_true_byHits_origin;
  vector<int>             *reco_daughter_PFP_true_byHits_parID;
  vector<int>             *reco_daughter_PFP_true_byHits_parPDG;
  vector<string>          *reco_daughter_PFP_true_byHits_process;
  vector<unsigned long>   *reco_daughter_PFP_true_byHits_sharedHits;
  vector<unsigned long>   *reco_daughter_PFP_true_byHits_emHits;
  vector<double>          *reco_daughter_PFP_true_byHits_len;
  vector<double>          *reco_daughter_PFP_true_byHits_startX;
  vector<double>          *reco_daughter_PFP_true_byHits_startY;
  vector<double>          *reco_daughter_PFP_true_byHits_startZ;
  vector<double>          *reco_daughter_PFP_true_byHits_endX;
  vector<double>          *reco_daughter_PFP_true_byHits_endY;
  vector<double>          *reco_daughter_PFP_true_byHits_endZ;
  vector<double>          *reco_daughter_PFP_true_byHits_startPx;
  vector<double>          *reco_daughter_PFP_true_byHits_startPy;
  vector<double>          *reco_daughter_PFP_true_byHits_startPz;
  vector<double>          *reco_daughter_PFP_true_byHits_startP;
  vector<double>          *reco_daughter_PFP_true_byHits_startE;
  vector<string>          *reco_daughter_PFP_true_byHits_endProcess;
  vector<double>          *reco_daughter_PFP_true_byHits_purity;
  vector<double>          *reco_daughter_PFP_true_byHits_completeness;
  vector<int>             *reco_daughter_PFP_true_byE_PDG;
  vector<double>          *reco_daughter_PFP_true_byE_len;
  vector<double>          *reco_daughter_PFP_true_byE_purity;
  vector<double>          *reco_daughter_PFP_true_byE_completeness; 
  // Beam Instrumentation Info
  Double_t        beam_inst_P;
  vector<double>  *beam_inst_TOF;
  vector<int>     *beam_inst_TOF_Chan;
  Double_t        beam_inst_X;
  Double_t        beam_inst_Y;
  Double_t        beam_inst_Z;
  Double_t        beam_inst_dirX;
  Double_t        beam_inst_dirY;
  Double_t        beam_inst_dirZ;
  vector<int>     *beam_inst_PDG_candidates;
  Int_t           beam_inst_nFibersP1;
  Int_t           beam_inst_nFibersP2;
  Int_t           beam_inst_nFibersP3;
  Int_t           beam_inst_nMomenta;
  Int_t           beam_inst_valid;
 
  // List of branches
  TBranch *b_run; //!
  TBranch *b_subrun; //!
  TBranch *b_event; //!
  TBranch *b_MC; //!
  TBranch *b_true_beam_PDG; //!
  TBranch *b_true_beam_mass; //!
  TBranch *b_true_beam_ID; //!
  TBranch *b_true_beam_endProcess; //!
  TBranch *b_true_beam_startX; //!
  TBranch *b_true_beam_startY; //!
  TBranch *b_true_beam_startZ; //!
  TBranch *b_true_beam_endX; //!
  TBranch *b_true_beam_endY; //!
  TBranch *b_true_beam_endZ; //!
  TBranch *b_true_beam_startP; //!
  TBranch *b_true_beam_startPx; //!
  TBranch *b_true_beam_startPy; //!
  TBranch *b_true_beam_startPz; //!
  TBranch *b_true_beam_endP; //!
  TBranch *b_true_beam_endPx; //!
  TBranch *b_true_beam_endPy; //!
  TBranch *b_true_beam_endPz; //!
  TBranch *b_true_beam_last_len; //!
  TBranch *b_true_beam_startDirX; //!
  TBranch *b_true_beam_startDirY; //!
  TBranch *b_true_beam_startDirZ; //!
  TBranch *b_true_beam_nElasticScatters; //!
  TBranch *b_true_beam_elastic_costheta; //!
  TBranch *b_true_beam_elastic_X; //!
  TBranch *b_true_beam_elastic_Y; //!
  TBranch *b_true_beam_elastic_Z; //!
  TBranch *b_true_beam_elastic_deltaE; //!
  TBranch *b_true_beam_elastic_IDE_edep; //!
  TBranch *b_true_beam_IDE_totalDep; //!
  TBranch *b_true_daughter_nPi0; //!
  TBranch *b_true_daughter_nPiPlus; //!
  TBranch *b_true_daughter_nProton; //!
  TBranch *b_true_daughter_nNeutron; //!
  TBranch *b_true_daughter_nPiMinus; //!
  TBranch *b_true_daughter_nNucleus; //!
  TBranch *b_true_beam_processes; //!
  TBranch *b_true_beam_incidentEnergies; //!
  TBranch *b_true_beam_interactingEnergy; //!
  TBranch *b_true_beam_slices; //!
  TBranch *b_true_beam_slices_deltaE; //!
  TBranch *b_true_beam_traj_X; //!
  TBranch *b_true_beam_traj_Y; //!
  TBranch *b_true_beam_traj_Z; //!
  TBranch *b_true_beam_traj_KE; //!
  TBranch *b_reco_beam_type; //!
  TBranch *b_reco_beam_trackID; //!
  TBranch *b_reco_beam_startX; //!
  TBranch *b_reco_beam_startY; //!
  TBranch *b_reco_beam_startZ; //!
  TBranch *b_reco_beam_endX; //!
  TBranch *b_reco_beam_endY; //!
  TBranch *b_reco_beam_endZ; //!
  TBranch *b_reco_beam_trackDirX; //!
  TBranch *b_reco_beam_trackDirY; //!
  TBranch *b_reco_beam_trackDirZ; //!
  TBranch *b_reco_beam_trackEndDirX; //!
  TBranch *b_reco_beam_trackEndDirY; //!
  TBranch *b_reco_beam_trackEndDirZ; //!
  TBranch *b_reco_beam_flipped; //!
  TBranch *b_reco_beam_len; //!
  TBranch *b_reco_beam_alt_len; //!
  TBranch *b_reco_beam_momByRange_proton; //!
  TBranch *b_reco_beam_momByRange_muon; //!
  TBranch *b_reco_beam_momByRange_alt_proton; //!
  TBranch *b_reco_beam_momByRange_alt_muon; //!
  TBranch *b_reco_beam_calo_startX; //!
  TBranch *b_reco_beam_calo_startY; //!
  TBranch *b_reco_beam_calo_startZ; //!
  TBranch *b_reco_beam_calo_endX; //!
  TBranch *b_reco_beam_calo_endY; //!
  TBranch *b_reco_beam_calo_endZ; //!
  TBranch *b_reco_beam_calo_startDirX; //!
  TBranch *b_reco_beam_calo_startDirY; //!
  TBranch *b_reco_beam_calo_startDirZ; //!
  TBranch *b_reco_beam_calo_endDirX; //!
  TBranch *b_reco_beam_calo_endDirY; //!
  TBranch *b_reco_beam_calo_endDirZ; //!
  TBranch *b_reco_beam_vertex_michel_score; //!
  TBranch *b_reco_beam_vertex_nHits; //!
  TBranch *b_reco_beam_dQdX_SCE; //!
  TBranch *b_reco_beam_dQ; //!
  TBranch *b_reco_beam_dEdX_SCE; //!
  TBranch *b_reco_beam_calibrated_dEdX_SCE; //!
  TBranch *b_reco_beam_resRange_SCE; //!
  TBranch *b_reco_beam_TrkPitch_SCE; //!
  TBranch *b_reco_beam_dQdX_NoSCE; //!
  TBranch *b_reco_beam_dEdX_NoSCE; //!
  TBranch *b_reco_beam_calibrated_dEdX_NoSCE; //!
  TBranch *b_reco_beam_resRange_NoSCE; //!
  TBranch *b_reco_beam_TrkPitch_NoSCE; //!
  TBranch *b_reco_beam_calo_wire; //!
  TBranch *b_reco_beam_calo_wire_z; //!
  TBranch *b_reco_beam_calo_tick; //!
  TBranch *b_reco_beam_calo_TPC; //!
  TBranch *b_reco_beam_passes_beam_cuts; //!
  TBranch *b_reco_beam_Chi2_proton; //!
  TBranch *b_reco_beam_Chi2_ndof; //!
  TBranch *b_reco_beam_PFP_ID; //!
  TBranch *b_reco_beam_PFP_nHits; //!
  TBranch *b_reco_beam_PFP_trackScore; //!
  TBranch *b_reco_beam_PFP_emScore; //!
  TBranch *b_reco_beam_PFP_michelScore; //!
  TBranch *b_reco_beam_PFP_trackScore_collection; //!
  TBranch *b_reco_beam_PFP_emScore_collection; //!
  TBranch *b_reco_beam_PFP_michelScore_collection; //!
  TBranch *b_reco_beam_allTrack_ID; //!
  TBranch *b_reco_beam_allTrack_beam_cuts; //!
  TBranch *b_reco_beam_allTrack_flipped; //!
  TBranch *b_reco_beam_allTrack_len; //!
  TBranch *b_reco_beam_allTrack_startX; //!
  TBranch *b_reco_beam_allTrack_startY; //!
  TBranch *b_reco_beam_allTrack_startZ; //!
  TBranch *b_reco_beam_allTrack_endX; //!
  TBranch *b_reco_beam_allTrack_endY; //!
  TBranch *b_reco_beam_allTrack_endZ; //!
  TBranch *b_reco_beam_allTrack_trackDirX; //!
  TBranch *b_reco_beam_allTrack_trackDirY; //!
  TBranch *b_reco_beam_allTrack_trackDirZ; //!
  TBranch *b_reco_beam_allTrack_trackEndDirX; //!
  TBranch *b_reco_beam_allTrack_trackEndDirY; //!
  TBranch *b_reco_beam_allTrack_trackEndDirZ; //!
  TBranch *b_reco_beam_allTrack_resRange; //!
  TBranch *b_reco_beam_allTrack_calibrated_dEdX; //!
  TBranch *b_reco_beam_allTrack_Chi2_proton; //!
  TBranch *b_reco_beam_allTrack_Chi2_ndof; //!
  TBranch *b_reco_beam_true_byHits_endProcess; //!
  TBranch *b_reco_beam_true_byHits_process; //!
  TBranch *b_reco_beam_true_byHits_origin; //!
  TBranch *b_reco_beam_true_byHits_PDG; //!
  TBranch *b_reco_beam_true_byHits_ID; //!
  TBranch *b_reco_beam_true_byHits_matched; //!
  TBranch *b_reco_beam_true_byHits_purity; //!
  TBranch *b_reco_beam_true_byHits_startE; //!
  TBranch *b_reco_beam_true_byHits_startP; //!
  TBranch *b_reco_beam_true_byHits_startPx; //!
  TBranch *b_reco_beam_true_byHits_startPy; //!
  TBranch *b_reco_beam_true_byHits_startPz; //!
  TBranch *b_reco_beam_true_byHits_endE; //!
  TBranch *b_reco_beam_true_byHits_endP; //!
  TBranch *b_reco_beam_true_byHits_endPx; //!
  TBranch *b_reco_beam_true_byHits_endPy; //!
  TBranch *b_reco_beam_true_byHits_endPz; //!
  TBranch *b_reco_beam_true_byE_startE; //!
  TBranch *b_reco_beam_true_byE_startP; //!
  TBranch *b_reco_beam_true_byE_startPx; //!
  TBranch *b_reco_beam_true_byE_startPy; //!
  TBranch *b_reco_beam_true_byE_startPz; //!
  TBranch *b_reco_beam_true_byE_endE; //!
  TBranch *b_reco_beam_true_byE_endP; //!
  TBranch *b_reco_beam_true_byE_endPx; //!
  TBranch *b_reco_beam_true_byE_endPy; //!
  TBranch *b_reco_beam_true_byE_endPz; //!
  TBranch *b_reco_beam_true_byE_endProcess; //!
  TBranch *b_reco_beam_true_byE_process; //!
  TBranch *b_reco_beam_true_byE_origin; //!
  TBranch *b_reco_beam_true_byE_PDG; //!
  TBranch *b_reco_beam_true_byE_ID; //!
  TBranch *b_reco_beam_true_byE_matched; //!
  TBranch *b_reco_beam_incidentEnergies; //!
  TBranch *b_reco_beam_interactingEnergy; //!
  TBranch *b_reco_track_startX; //!
  TBranch *b_reco_track_startY; //!
  TBranch *b_reco_track_startZ; //!
  TBranch *b_reco_track_endX; //!
  TBranch *b_reco_track_endY; //!
  TBranch *b_reco_track_endZ; //!
  TBranch *b_reco_track_michel_score; //!
  TBranch *b_reco_track_nHits; //!
  TBranch *b_reco_track_ID; //!
  TBranch *b_true_beam_daughter_ID; //!
  TBranch *b_true_beam_daughter_PDG; //!
  TBranch *b_true_beam_daughter_len; //!
  TBranch *b_true_beam_daughter_startX; //!
  TBranch *b_true_beam_daughter_startY; //!
  TBranch *b_true_beam_daughter_startZ; //!
  TBranch *b_true_beam_daughter_endX; //!
  TBranch *b_true_beam_daughter_endY; //!
  TBranch *b_true_beam_daughter_endZ; //!
  TBranch *b_true_beam_daughter_startP; //!
  TBranch *b_true_beam_daughter_startPx; //!
  TBranch *b_true_beam_daughter_startPy; //!
  TBranch *b_true_beam_daughter_startPz; //!
  TBranch *b_true_beam_daughter_Process; //!
  TBranch *b_true_beam_daughter_endProcess; //!
  TBranch *b_true_beam_grand_daughter_ID; //!
  TBranch *b_true_beam_grand_daughter_parID; //!
  TBranch *b_true_beam_grand_daughter_PDG; //!
  TBranch *b_true_beam_grand_daughter_Process; //!
  TBranch *b_true_beam_grand_daughter_endProcess; //!
  TBranch *b_true_beam_Pi0_decay_ID; //!
  TBranch *b_true_beam_Pi0_decay_parID; //!
  TBranch *b_true_beam_Pi0_decay_PDG; //!
  TBranch *b_true_beam_Pi0_decay_startX; //!
  TBranch *b_true_beam_Pi0_decay_startY; //!
  TBranch *b_true_beam_Pi0_decay_startZ; //!
  TBranch *b_true_beam_Pi0_decay_startP; //!
  TBranch *b_true_beam_Pi0_decay_startPx; //!
  TBranch *b_true_beam_Pi0_decay_startPy; //!
  TBranch *b_true_beam_Pi0_decay_startPz; //!
  TBranch *b_true_beam_Pi0_decay_len; //!
  TBranch *b_reco_daughter_PFP_ID; //!
  TBranch *b_reco_daughter_PFP_nHits; //!
  TBranch *b_reco_daughter_PFP_trackScore; //!
  TBranch *b_reco_daughter_PFP_emScore; //!
  TBranch *b_reco_daughter_PFP_michelScore; //!
  TBranch *b_reco_daughter_PFP_nHits_collection; //!
  TBranch *b_reco_daughter_PFP_trackScore_collection; //!
  TBranch *b_reco_daughter_PFP_emScore_collection; //!
  TBranch *b_reco_daughter_PFP_michelScore_collection; //!
  TBranch *b_reco_daughter_allTrack_ID; //!
  TBranch *b_reco_daughter_allTrack_Theta; //!
  TBranch *b_reco_daughter_allTrack_Phi; //!
  TBranch *b_reco_daughter_allTrack_len; //!
  TBranch *b_reco_daughter_allTrack_alt_len; //!
  TBranch *b_reco_daughter_allTrack_momByRange_proton; //!
  TBranch *b_reco_daughter_allTrack_momByRange_muon; //!
  TBranch *b_reco_daughter_allTrack_momByRange_alt_proton; //!
  TBranch *b_reco_daughter_allTrack_momByRange_alt_muon; //!
  TBranch *b_reco_daughter_allTrack_startX; //!
  TBranch *b_reco_daughter_allTrack_startY; //!
  TBranch *b_reco_daughter_allTrack_startZ; //!
  TBranch *b_reco_daughter_allTrack_endX; //!
  TBranch *b_reco_daughter_allTrack_endY; //!
  TBranch *b_reco_daughter_allTrack_endZ; //!
  TBranch *b_reco_daughter_allTrack_vertex_michel_score; //!
  TBranch *b_reco_daughter_allTrack_vertex_nHits; //!
  TBranch *b_reco_daughter_allTrack_dQdX_SCE; //!
  TBranch *b_reco_daughter_allTrack_dEdX_SCE; //!
  TBranch *b_reco_daughter_allTrack_calibrated_dEdX_SCE; //!
  TBranch *b_reco_daughter_allTrack_resRange_SCE; //!
  TBranch *b_reco_daughter_allTrack_Chi2_proton; //!
  TBranch *b_reco_daughter_allTrack_Chi2_ndof; //!
  TBranch *b_reco_daughter_allTrack_calibrated_dEdX_SCE_plane0; //!
  TBranch *b_reco_daughter_allTrack_resRange_plane0; //!
  TBranch *b_reco_daughter_allTrack_Chi2_proton_plane0; //!
  TBranch *b_reco_daughter_allTrack_Chi2_ndof_plane0; //!
  TBranch *b_reco_daughter_allTrack_calibrated_dEdX_SCE_plane1; //!
  TBranch *b_reco_daughter_allTrack_resRange_plane1; //!
  TBranch *b_reco_daughter_allTrack_Chi2_proton_plane1; //!
  TBranch *b_reco_daughter_allTrack_Chi2_ndof_plane1; //!
  TBranch *b_reco_daughter_allShower_ID; //!
  TBranch *b_reco_daughter_allShower_len; //!
  TBranch *b_reco_daughter_allShower_startX; //!
  TBranch *b_reco_daughter_allShower_startY; //!
  TBranch *b_reco_daughter_allShower_startZ; //!
  TBranch *b_reco_daughter_allShower_dirX; //!
  TBranch *b_reco_daughter_allShower_dirY; //!
  TBranch *b_reco_daughter_allShower_dirZ; //!
  TBranch *b_reco_daughter_allShower_energy; //!
  TBranch *b_reco_daughter_PFP_true_byHits_PDG; //!
  TBranch *b_reco_daughter_PFP_true_byHits_ID; //!
  TBranch *b_reco_daughter_PFP_true_byHits_origin; //!
  TBranch *b_reco_daughter_PFP_true_byHits_parID; //!
  TBranch *b_reco_daughter_PFP_true_byHits_parPDG; //!
  TBranch *b_reco_daughter_PFP_true_byHits_process; //!
  TBranch *b_reco_daughter_PFP_true_byHits_sharedHits; //!
  TBranch *b_reco_daughter_PFP_true_byHits_emHits; //!
  TBranch *b_reco_daughter_PFP_true_byHits_len; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startX; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startY; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startZ; //!
  TBranch *b_reco_daughter_PFP_true_byHits_endX; //!
  TBranch *b_reco_daughter_PFP_true_byHits_endY; //!
  TBranch *b_reco_daughter_PFP_true_byHits_endZ; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startPx; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startPy; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startPz; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startP; //!
  TBranch *b_reco_daughter_PFP_true_byHits_startE; //!
  TBranch *b_reco_daughter_PFP_true_byHits_endProcess; //!
  TBranch *b_reco_daughter_PFP_true_byHits_purity; //!
  TBranch *b_reco_daughter_PFP_true_byHits_completeness; //!
  TBranch *b_reco_daughter_PFP_true_byE_PDG; //!
  TBranch *b_reco_daughter_PFP_true_byE_len; //!
  TBranch *b_reco_daughter_PFP_true_byE_purity; //!
  TBranch *b_reco_daughter_PFP_true_byE_completeness; //! 
  TBranch *b_beam_inst_P; //!
  TBranch *b_beam_inst_TOF; //!
  TBranch *b_beam_inst_TOF_Chan; //!
  TBranch *b_beam_inst_X; //!
  TBranch *b_beam_inst_Y; //!
  TBranch *b_beam_inst_Z; //!
  TBranch *b_beam_inst_dirX; //!
  TBranch *b_beam_inst_dirY; //!
  TBranch *b_beam_inst_dirZ; //!
  TBranch *b_beam_inst_PDG_candidates; //!
  TBranch *b_beam_inst_nFibersP1; //!
  TBranch *b_beam_inst_nFibersP2; //!
  TBranch *b_beam_inst_nFibersP3; //!
  TBranch *b_beam_inst_nMomenta; //!
  TBranch *b_beam_inst_valid; //!                                                                             
  
  
  // Header's
  Int_t EventTime; 
  Int_t TriggerWord; 
  Float_t POTPerSpill; 
}; 



#endif

