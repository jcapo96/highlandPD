/*
   converter for michelremovingtrees 

   M. García December 2022
*/

#ifndef MichelRemovingTreeConverter_h
#define MichelRemovingTreeConverter_h

#include "pdBaseConverter.hxx"

class MichelRemovingTreeConverter: public pdBaseConverter{
  
public:
  
  MichelRemovingTreeConverter(const std::string& name);
  virtual ~MichelRemovingTreeConverter(){}
  
  virtual void InitializeVariables();
  virtual void SetBranchAddresses();

  // ----------------------------
  
  virtual void FillEventInfo(AnaEventInfo* info);  
  virtual void FillTrueInfo(AnaSpill* spill){return;} //no true info needed so far
  virtual void FillBeamInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBeamPD* beam){return;} //no beam info needed so far
  virtual void FillBunchInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBunch* bunch, AnaBeamPD* beam);

  void FillParticleInfo(AnaParticlePD* part, const int itrk);
  bool IsUsableHit(const int itrk, const int ihit);
  

protected:
  
  const static int _MaxTracks = 30;
  const static int _NPlanes = 3;
  const static int _MaxHits = 3000;


  //Declaration of leaf types
  //MetaData
  bool     isData;
  Int_t    run;                  
  Int_t    subrun;               
  Int_t    event;
  // Double_t evttime; 
  // Int_t    year_month_date;
  // Int_t    hour_min_sec;

  //Number of tracks
  Int_t    cross_trks; //are the total number of good T0 tagged tracks stopping + crossing  
  // Int_t    stopping_trks;
  // Int_t    all_trks;
  // Int_t    unbroken_trks; //these are unbroken stopping tracks

  //Tracks information
  Float_t  trackthetaxz[_MaxTracks];
  Float_t  trackthetayz[_MaxTracks];
  Float_t  trkstartx[_MaxTracks];
  Float_t  trkstarty[_MaxTracks];
  Float_t  trkstartz[_MaxTracks];
  Float_t  trkendx[_MaxTracks];
  Float_t  trkendy[_MaxTracks];
  Float_t  trkendz[_MaxTracks];
  Float_t  trklen[_MaxTracks];
  Int_t    TrkID[_MaxTracks]; 
  // Float_t  trkstartcosxyz[_MaxTracks][_NPlanes];
  // Float_t  trkendcosxyz[_MaxTracks][_NPlanes];
  Int_t    ntrkhits[_MaxTracks][_NPlanes];
  Float_t  trkdqdx[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkdedx[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkresrange[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkhitx[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkhity[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkhitz[_MaxTracks][_NPlanes][_MaxHits];
  Float_t  trkpitch[_MaxTracks][_NPlanes][_MaxHits];
  // Float_t  peakT_max[_MaxTracks];
  // Float_t  peakT_min[_MaxTracks];
  // Float_t  dist_min[_MaxTracks];
  // Int_t    adjacent_hits[_MaxTracks];
  // Int_t    lastwire[_MaxTracks];
  // Int_t    endtpc[_MaxTracks];
  // Float_t  lastpeakt[_MaxTracks];

  //Tracks true info → not used
  // float    true_trkstartx[_MaxTracks];
  // float    true_trkstarty[_MaxTracks];
  // float    true_trkstartz[_MaxTracks];
  // float    true_trkendx[_MaxTracks];
  // float    true_trkendy[_MaxTracks];
  // float    true_trkendz[_MaxTracks];
 
  //List of branches
  TBranch *b_isData; //!
  TBranch *b_run; //!                  
  TBranch *b_subrun; //!               
  TBranch *b_event; //!
  // TBranch *b_evttime; //! 
  // TBranch *b_year_month_date; //!
  // TBranch *b_hour_min_sec; //!
  TBranch *b_cross_trks; //! //are the total number of good T0 tagged tracks stopping + crossing  
  // TBranch *b_stopping_trks; //!
  // TBranch *b_all_trks; //!
  // TBranch *b_unbroken_trks; //! //these are unbroken stopping tracks
  TBranch *b_trackthetaxz; //!
  TBranch *b_trackthetayz; //!
  TBranch *b_trkstartx; //!
  TBranch *b_trkstarty; //!
  TBranch *b_trkstartz; //!
  TBranch *b_trkendx; //!
  TBranch *b_trkendy; //!
  TBranch *b_trkendz; //!
  TBranch *b_trklen; //!
  TBranch *b_TrkID; //! 
  // TBranch *b_trkstartcosxyz; //!
  //TBranch *b_trkendcosxyz; //!
  TBranch *b_ntrkhits; //!
  TBranch *b_trkdqdx; //!
  TBranch *b_trkdedx; //!
  TBranch *b_trkresrange; //!
  TBranch *b_trkhitx; //!
  TBranch *b_trkhity; //!
  TBranch *b_trkhitz; //!
  TBranch *b_trkpitch; //!
  // TBranch *b_peakT_max; //!
  // TBranch *b_peakT_min; //!
  // TBranch *b_dist_min; //!
  // TBranch *b_adjacent_hits; //!
  // TBranch *b_lastwire; //!
  // TBranch *b_endtpc; //!
  // TBranch *b_lastpeakt; //!
  // TBranch *b_true_trkstartx; //!
  // TBranch *b_true_trkstarty; //!
  // TBranch *b_true_trkstartz; //!
  // TBranch *b_true_trkendx; //!
  // TBranch *b_true_trkendy; //!
  // TBranch *b_true_trkendz; //!
};  



#endif

