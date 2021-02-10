#ifndef CNNUtils_h
#define CNNUtils_h


#include "PointIdAlg.hxx"
#include "pdDataClasses.hxx"

class CNNUtils{
public:
  
  
  CNNUtils();
  
  
  void produce(std::vector<AnaHitPD>& hits); //override;
  void FillHits(const AnaEventPD& evt, std::vector<AnaHitPD*>& hitPtrList);

  
private:    
  size_t fBatchSize;
  nnet::PointIdAlg fPointIdAlg;
};
 
  
#endif
