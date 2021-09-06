#ifndef pdMiniTreeConverter_h
#define pdMiniTreeConverter_h

#include "pdDataClasses.hxx"
#include "HighlandMiniTreeConverter.hxx"

class pdMiniTreeConverter: public HighlandMiniTreeConverter{

 public:
  pdMiniTreeConverter(const std::string& tree_path="MiniTree", bool readRooTrackerVtx=false):HighlandMiniTreeConverter(tree_path, readRooTrackerVtx){}
  virtual ~pdMiniTreeConverter(){}

  // Create the appropriate spill instance
  virtual AnaSpillB*         MakeSpill() { return new AnaSpill(); }
  //  virtual AnaEventB*         MakeEvent() { return new AnaEventPD(); }


  Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill);
  
};

#endif


