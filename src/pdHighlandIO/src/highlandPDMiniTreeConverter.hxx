#ifndef highlandPDMiniTreeConverter_h
#define highlandPDMiniTreeConverter_h

#include "pdDataClasses.hxx"
#include "HighlandMiniTreeConverter.hxx"

/// Creates the appropriate AnaSpillB type. The rest of the work is done by the base converter

class highlandPDMiniTreeConverter: public MiniTreeConverter{

 public:
  highlandPDMiniTreeConverter(const std::string& tree_path="MiniTree", bool readRooTrackerVtx=false):MiniTreeConverter(tree_path, readRooTrackerVtx){}
  virtual ~highlandPDMiniTreeConverter(){}

  // Create the appropriate spill instance
  virtual AnaSpillB*         MakeSpill() { return new AnaSpillPD(); }
  virtual AnaEventB*         MakeEvent() { return new AnaEventPD(); }


//*****************************************************************************
Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill){
//*****************************************************************************

  // Read contents of entry.
  if (!fChain) return 0;
  
  // Create the appropriate instance of Spill
  spill = MakeSpill();
  
  // cast to AnaSpillB
  AnaSpillPD* _spill2 = static_cast<AnaSpillPD*>(spill);

  std::cout <<"anselmo 1" << std::endl;
  // Set the branch address
  fChain->SetBranchAddress("Spill", &spill);
  std::cout <<"anselmo 2" << std::endl;
  // Print the current file
  static std::string currentfilename=""; 
  std::string filename =minitree->GetFile()->GetName();   
  if( filename != currentfilename ) {
    std::cout << " Running on file: " << filename << std::endl; 
    currentfilename = filename;
  }

  // get a new entry from the flat tree. entry_temp >0 when succesfull
  Int_t entry_temp = minitree->GetEntry(entry);
  std::cout <<"anselmo 3" << std::endl;
  return 0;
  if (_spill->EventInfo){
    if (_readRooTrackerVtx){
      bool sIsMC    = _spill->EventInfo->IsMC;
      Int_t sEvt    = _spill->EventInfo->Event;
      Int_t sSubrun = _spill->EventInfo->SubRun;
      Int_t sRun    = _spill->EventInfo->Run;
      
      // sEvt should be positive since sEvt=-999 is used for the last flatree entry
      if (entry_temp>0 && sIsMC && (fGenie || fNeut) && sEvt>=0) {
        // Loop over RooTrackerVtx entries until we get the same run, subrun and event numbers
        // In general we will get the right entry in the first iteration, but just in case
        do{       
          if      (NRooTrackerVTX)  NRooTrackerVTX->GetEntry(_entry_roo);
          else if (GRooTrackerVTX)  GRooTrackerVTX->GetEntry(_entry_roo);
          if ((RunID> sRun                                      ) ||
              (RunID==sRun && SubrunID> sSubrun                 ) ||
              (RunID==sRun && SubrunID==sSubrun && EventID>sEvt ))
            _entry_roo--;  
          else
            _entry_roo++;  
        }while(EventID!=sEvt || RunID!=sRun || SubrunID!=sSubrun);
      }
    }
    
    // Copy vectors into arrays
    _spill->CopyVectorsIntoArrays();
    // Redo reco-reco and reco-truth links (many not pressent in MiniTree)
    _spill->RedoLinks();
  }
  else
    entry_temp=0;


  // Increment entry number
  entry++;
  if (entry%10000==0 || entry == _nentries || (entry%1000==0 && entry<10000) )
    std::cout << "entry: " << entry << " of " << _nentries << " (" << (100*entry/_nentries) << "%)" << std::endl;

  // Load the geometry for this spill(999 is the default value in BaseDataClasses). TODO DUNE
  //  if (_spill->GeomID!=999)
    //    ND::hgman().LoadGeometry(filename,(Int_t)_spill->GeomID,"geom");

  return entry_temp;
}


};

#endif


