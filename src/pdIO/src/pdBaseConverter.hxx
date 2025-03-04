/*
   general base converter for ProtoDUNE-SP flat trees

   A. Cervera September 2021
*/

#ifndef pdBaseConverter_h
#define pdBaseConverter_h

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
#include "pdDataClasses.hxx"

using namespace std;

class pdBaseConverter: public InputConverter{

 public:

  pdBaseConverter(const std::string& name);
  virtual ~pdBaseConverter();

  virtual bool Initialize();
  virtual Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill);
  Int_t GetEvent(Long64_t& entry, AnaEventC*& event){(void)entry;(void)event; return 0;}

  /// Record the POT for the current spill, based on information in the AnaBeam
  /// member of the current AnaSpill.
  void IncrementPOTBySpill(){return;}

  virtual Int_t ReadEntry(Long64_t& entry);
  virtual bool AddFileToTChain(const std::string& inputString);
  virtual void InitializeVariables()=0;
  virtual void SetBranchAddresses()=0;
  
  //----------------
  virtual AnaSpill*           MakeSpill()       { return new AnaSpill(); }
  virtual AnaBunch*           MakeBunch()       { return new AnaBunch(); }
  virtual AnaBeamPD*          MakeBeam()        { return new AnaBeamPD(); }
  virtual AnaDataQuality*     MakeDataQuality() { return new AnaDataQuality(); }
  virtual AnaEventInfoPD*     MakeEventInfo()   { return new AnaEventInfoPD(); }
  virtual AnaTrigger*         MakeTrigger()     { return new AnaTrigger(); }

  virtual AnaTrueParticlePD*  MakeTrueParticle(){ return new AnaTrueParticlePD(); }
  virtual AnaTrueVertex*      MakeTrueVertex()  { return new AnaTrueVertex(); }
  virtual AnaParticlePD*      MakeParticle()    { return new AnaParticlePD(); }
  // ----------------------------

  virtual void FillInfo(AnaSpill* spill);
  virtual void FillDQInfo(AnaDataQuality* dq);
  
  // ---- pure virtual methods to be inplemented in derived class  ----
  virtual void FillTrueInfo(AnaSpill* spill) = 0;
  virtual void FillBeamInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBeamPD* beam) = 0;
  virtual void FillBunchInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBunch* bunch, AnaBeamPD* beam) = 0;
  virtual void FillEventInfo(AnaEventInfoPD* info) = 0;
  //-----------------

protected:

  AnaSpill* _spill;
  
  std::string _previousFile;
  Int_t _previousRunID;
  Int_t _previousSubrunID;
  Int_t _previousRefEventID;

 protected:

  Bool_t _isMC;
  std::string _softwareVersion;

  Bool_t _FillEventInfo;
  Bool_t _FillDQInfo;
  Bool_t _FillBeamInfo;
  Bool_t _FillBunchInfo;
  Bool_t _FillTrueInfo;
}; 

#endif

