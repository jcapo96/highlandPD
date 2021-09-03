#ifndef highlandPDMiniTreeConverter_h
#define highlandPDMiniTreeConverter_h

#include "pdDataClasses.hxx"
#include "HighlandMiniTreeConverter.hxx"

class highlandPDMiniTreeConverter: public HighlandMiniTreeConverter{

 public:
  highlandPDMiniTreeConverter(const std::string& tree_path="MiniTree", bool readRooTrackerVtx=false):HighlandMiniTreeConverter(tree_path, readRooTrackerVtx){}
  virtual ~highlandPDMiniTreeConverter(){}

  // Create the appropriate spill instance
  virtual AnaSpillB*         MakeSpill() { return new AnaSpill(); }
  //  virtual AnaEventB*         MakeEvent() { return new AnaEventPD(); }


  Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill);
  
};

#endif


