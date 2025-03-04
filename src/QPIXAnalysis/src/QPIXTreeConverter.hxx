#ifndef QPIXTreeConverter_h
#define QPIXTreeConverter_h

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
#include "QPIXDataClasses.hxx"

using namespace std;

class QPIXTreeConverter: public InputConverter{

 public:

  QPIXTreeConverter(const std::string& name);
  virtual ~QPIXTreeConverter();

  virtual bool Initialize();
  virtual Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill);
  Int_t GetEvent(Long64_t& entry, AnaEventC*& event){(void)entry;(void)event; return 0;}

  /// Record the POT for the current spill, based on information in the AnaBeam
  /// member of the current AnaSpill.
  void IncrementPOTBySpill(){return;}

  virtual Int_t ReadEntry(Long64_t& entry);
  virtual bool AddFileToTChain(const std::string& inputString);
  virtual void InitializeVariables();
  virtual void SetBranchAddresses();
  
  //----------------
  virtual AnaSpill*           MakeSpill()       { return new AnaSpill(); }
  virtual AnaBunch*           MakeBunch()       { return new AnaBunch(); }
  virtual AnaEventInfoPD*     MakeEventInfo()   { return new AnaEventInfoQPIX(); }
  virtual AnaBeamPD*          MakeBeam()        { return new AnaBeamPD(); }
  virtual AnaDataQuality*     MakeDataQuality() { return new AnaDataQuality(); }

  virtual AnaTrueParticlePD*  MakeTrueParticle(){ return new AnaTrueParticlePD(); }
  virtual AnaTrueVertex*      MakeTrueVertex()  { return new AnaTrueVertex(); }
  // ----------------------------

  virtual void FillInfo(AnaSpill* spill);
  
  virtual void FillTrueInfo(AnaSpill* spill);
  virtual void FillEventInfo(AnaEventInfoQPIX* info);
  //-----------------

  void GetSolarChain(AnaEventInfoQPIX* info);
  void GetNuReaction(AnaEventInfoQPIX* info);
  
protected:

  std::string _filename;

  AnaSpill* _spill;

  // Declaration of leaf types
  // Meta data
  int _qpix_run;
  int _qpix_event;
  int _number_particles;
  int _is_background;
  int _has_potasium;
  int _has_gamma;   

  //all particles of the event
  std::vector<int>    *_particle_track_id;
  std::vector<int>    *_particle_parent_track_id;
  std::vector<int>    *_particle_pdg_code;
  std::vector<double> *_particle_initial_energy;
  std::vector<double> *_particle_mass;
  std::vector<double> *_particle_charge;
  std::vector<int>    *_particle_process_key;
  std::vector<double> *_particle_initial_x;
  std::vector<double> *_particle_initial_y;
  std::vector<double> *_particle_initial_z;
  std::vector<double> *_particle_initial_t;
  std::vector<double> *_particle_initial_px;
  std::vector<double> *_particle_initial_py;
  std::vector<double> *_particle_initial_pz;
  std::vector<int>    *_particle_number_daughters;
  
  // all particle hits
  std::vector<int>    *_hit_track_id;
  std::vector<double> *_hit_start_x;
  std::vector<double> *_hit_start_y;
  std::vector<double> *_hit_start_z;
  std::vector<double> *_hit_start_t;
  std::vector<double> *_hit_end_x;
  std::vector<double> *_hit_end_y;
  std::vector<double> *_hit_end_z;
  std::vector<double> *_hit_end_t;
  std::vector<double> *_hit_length;
  std::vector<double> *_hit_energy_deposit;
  std::vector<int >   *_hit_process_key;

  // vectors for initial generator particles
  std::vector<double> *_generator_initial_particle_energy;
  std::vector<int>    *_generator_initial_particle_pdg_code;
  std::vector<double> *_generator_initial_particle_mass;
  std::vector<double> *_generator_initial_particle_charge;
  std::vector<double> *_generator_initial_particle_x;
  std::vector<double> *_generator_initial_particle_y;
  std::vector<double> *_generator_initial_particle_z;
  std::vector<double> *_generator_initial_particle_t;
  std::vector<double> *_generator_initial_particle_px;
  std::vector<double> *_generator_initial_particle_py;
  std::vector<double> *_generator_initial_particle_pz;

  // vectors for final generator particles
  std::vector<double> *_generator_final_particle_energy;
  std::vector<int>    *_generator_final_particle_pdg_code;
  std::vector<double> *_generator_final_particle_mass;
  std::vector<double> *_generator_final_particle_charge;
  std::vector<double> *_generator_final_particle_x;
  std::vector<double> *_generator_final_particle_y;
  std::vector<double> *_generator_final_particle_z;
  std::vector<double> *_generator_final_particle_t;
  std::vector<double> *_generator_final_particle_px;
  std::vector<double> *_generator_final_particle_py;
  std::vector<double> *_generator_final_particle_pz;           

  // waveforms and detected photons
  std::vector<int> *_wf[7];
  std::vector<int> *_detected_photons;

  //TBranches
  TBranch* b_qpix_run;//!
  TBranch* b_qpix_event;//!
  TBranch* b_number_particles;//!
  TBranch* b_is_background;//!
  TBranch* b_has_potasium;//!
  TBranch* b_has_gamma;//!
  TBranch* b_particle_track_id;//!
  TBranch* b_particle_parent_track_id;//!
  TBranch* b_particle_pdg_code;//!
  TBranch* b_particle_initial_energy;//!
  TBranch* b_particle_mass;//!
  TBranch* b_particle_charge;//!
  TBranch* b_particle_process_key;//!
  TBranch* b_particle_initial_x;//!
  TBranch* b_particle_initial_y;//!
  TBranch* b_particle_initial_z;//!
  TBranch* b_particle_initial_t;//!
  TBranch* b_particle_initial_px;//!
  TBranch* b_particle_initial_py;//!
  TBranch* b_particle_initial_pz;//!
  TBranch* b_particle_number_daughters;//!
  TBranch* b_hit_track_id;//!
  TBranch* b_hit_start_x;//!
  TBranch* b_hit_start_y;//!
  TBranch* b_hit_start_z;//!
  TBranch* b_hit_start_t;//!
  TBranch* b_hit_end_x;//!
  TBranch* b_hit_end_y;//!
  TBranch* b_hit_end_z;//!
  TBranch* b_hit_end_t;//!
  TBranch* b_hit_length;//!
  TBranch* b_hit_energy_deposit;//!
  TBranch* b_hit_process_key;//!
  TBranch* b_generator_initial_particle_energy;//!
  TBranch* b_generator_initial_particle_pdg_code;//!
  TBranch* b_generator_initial_particle_mass;//!
  TBranch* b_generator_initial_particle_charge;//!
  TBranch* b_generator_initial_particle_x;//!
  TBranch* b_generator_initial_particle_y;//!
  TBranch* b_generator_initial_particle_z;//!
  TBranch* b_generator_initial_particle_t;//!
  TBranch* b_generator_initial_particle_px;//!
  TBranch* b_generator_initial_particle_py;//!
  TBranch* b_generator_initial_particle_pz;//!
  TBranch* b_generator_final_particle_energy;//!
  TBranch* b_generator_final_particle_pdg_code;//!
  TBranch* b_generator_final_particle_mass;//!
  TBranch* b_generator_final_particle_charge;//!
  TBranch* b_generator_final_particle_x;//!
  TBranch* b_generator_final_particle_y;//!
  TBranch* b_generator_final_particle_z;//!
  TBranch* b_generator_final_particle_t;//!
  TBranch* b_generator_final_particle_px;//!
  TBranch* b_generator_final_particle_py;//!
  TBranch* b_generator_final_particle_pz;//!           
  TBranch* b_wf[7];//!
  TBranch* b_detected_photons;//!
}; 

#endif

