#include "pdBaseConverter.hxx"
#include "InputManager.hxx"
#include "BasicUtils.hxx"
#include "HighlandAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include "CalorimetryAlg.hxx"


//********************************************************************
pdBaseConverter::pdBaseConverter(const std::string& name):InputConverter(name){
//********************************************************************

  //constructor
  _spill = NULL;

  _isMC = false;
  _softwareVersion = "";

  _previousFile = "";
  _previousRunID = -1;
  _previousSubrunID = -1;
  _previousRefEventID = -1;
}

//********************************************************************
bool pdBaseConverter::Initialize(){
//********************************************************************


  AddChain(_treeName);
  eventsTree = GetChain(_treeName);
  
  fChain = eventsTree;
  
  // Set object pointer
  InitializeVariables();

  // Set branch addresses and branch pointers
  if (!fChain) return false;
  fCurrent = -1;
  SetBranchAddresses();
  
  return true;
}

//********************************************************************
pdBaseConverter::~pdBaseConverter(){
//********************************************************************
  
  //destructor 
  if(!fChain)return;

  if(eventsTree)delete eventsTree->GetCurrentFile();
}

//****************************************************************************
bool pdBaseConverter::AddFileToTChain(const std::string& inputString){
//****************************************************************************

  std::cout << "pdBaseConverter::AddFileToTChain(). Adding file: " << inputString << std::endl;

  // Chain only the directories we are interested in
  if(eventsTree)eventsTree->AddFile(inputString.c_str());
  
  // Read one entry from the tree tree such that Run and Subrun are available
  eventsTree->GetEntry(1);


  AnaEventInfo* evtInfo = MakeEventInfo();

  // general event info
  FillEventInfo(evtInfo);

  
  // Make sure the current file has not the same run and subrun number as the previous
  if (_previousRunID==evtInfo->Run &&  _previousSubrunID==evtInfo->SubRun && _previousRefEventID>= evtInfo->Event){
    std::cout << "-------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "pdBaseConverter::AddFileToTChain(). Current file has the same run and subrun as the previous" << std::endl;
    std::cout << "                                           and no higher event number !!!" << std::endl;
    std::cout << "   - this file:     " << inputString << std::endl;
    std::cout << "   - previous file: " << _previousFile << std::endl;
    std::cout << "Please verify the input file list !!!" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------" << std::endl;
    exit(1);
  }
  
  // The previous attributes
  _previousFile         = inputString;
  _previousRunID        = evtInfo->Run;
  _previousSubrunID     = evtInfo->SubRun;
  _previousRefEventID   = evtInfo->Event;
  
  // Set the data/MC mode and return false when mixing data and MC files
  _isMC = (bool)evtInfo->IsMC;
  if(!header().SetIsMC(_isMC)) return false;

  _softwareVersion = "v09_24"; //not sure about this but not important


  delete evtInfo;
  
  // Sets the software version for this file
  return header().SetSoftwareVersion(_softwareVersion);
}


//*****************************************************************************
Int_t pdBaseConverter::ReadEntries(Long64_t& entry) {
//*****************************************************************************
  
  Int_t entry_temp = eventsTree->GetEntry(entry);

  return entry_temp;
}

//*****************************************************************************
Int_t pdBaseConverter::GetSpill(Long64_t& entry, AnaSpillC*& spill){
//*****************************************************************************

  Int_t entry_temp = ReadEntries(entry);

  if(entry_temp > 0){
    
    // Create an instance of the Spill
    spill = MakeSpill();
    
    // Cast it to AnaSpill
    _spill = static_cast<AnaSpill*>(spill);
    
    // Fill the EventModel
    FillInfo(_spill);
  }
  else{
    std::cout << "Failed in reading entry " << entry << std::endl;
  }

  entry++;

  return entry_temp;
}

//*****************************************************************************
void pdBaseConverter::FillInfo(AnaSpill* spill){
//*****************************************************************************

  spill->EventInfo   = MakeEventInfo();
  spill->DataQuality = MakeDataQuality();
  spill->Beam        = MakeBeam();

  // general event info
  FillEventInfo(static_cast<AnaEventInfo*>(spill->EventInfo));
  
  // data quality info
  FillDQInfo(static_cast<AnaDataQuality*>(spill->DataQuality));

  // True information
  FillTrueInfo(spill);

  // beam related information (must be after true info)
  FillBeamInfo(spill->TrueParticles, static_cast<AnaBeamPD*>(spill->Beam));
  
  // All information about each bunch (only one in pd) (reco info)
  AnaBunch* bunch = MakeBunch();
  spill->Bunches.push_back(bunch);
  FillBunchInfo(spill->TrueParticles, bunch, static_cast<AnaBeamPD*>(spill->Beam));
}

//*****************************************************************************
void pdBaseConverter::FillDQInfo(AnaDataQuality* dq){
//*****************************************************************************

  dq->GoodDaq = true;
}


