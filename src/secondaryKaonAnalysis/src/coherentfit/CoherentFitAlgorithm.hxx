#ifndef CoherentFitAlgorithm_hxx
#define CoherentFitAlgorithm_hxx

#include "CoherentFit.hxx"

class CoherentFitAlgorithm{
public:
  CoherentFitAlgorithm(int argc, char *argv[]);
  virtual ~CoherentFitAlgorithm(){}

  void Execute();

  void ProcessMC();
  void ProcessData();

private:

  CoherentFit* fdata;
  CoherentFit* fmc;

  double fRMIN; 
  double fRMAX;
  double fSTEP;
  int    fNBINS;
  double fChi2Cut;
};

#endif
