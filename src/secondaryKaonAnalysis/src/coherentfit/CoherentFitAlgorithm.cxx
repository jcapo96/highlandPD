#include "CoherentFitAlgorithm.hxx"
#include <iostream>

//********************************************************************
CoherentFitAlgorithm::CoherentFitAlgorithm(int argc, char *argv[]){
//********************************************************************

  //primary initializations
  fRMIN    = 1; 
  fRMAX    = 61;
  fSTEP    = 2;
  fNBINS   = 30;
  fChi2Cut = 84;
  
  std::string dfile;
  std::string mcfile;

  for(;;){
    int c = getopt(argc, argv, "d:m:");
    if(c < 0)
      break;
    
    switch(c){
    case 'd':{
      dfile = optarg;
      break;
    }
    case 'm':{
      mcfile = optarg;
      break;
    }
    default: {
      std::cout << "MONIATO" << std::endl;
      std::exit(1);
    }
    }
  }

  if(dfile == "" || mcfile == ""){
    std::cout << "MONIATO" << std::endl;
    std::exit(1);
  }
  
  fdata = new CoherentFit(dfile);
  fmc   = new CoherentFit(mcfile);
}

//********************************************************************
void CoherentFitAlgorithm::Execute(){
//********************************************************************

  ProcessMC();
  ProcessData();
  fdata->GenerateHistograms(fRMIN,fRMAX,fSTEP,fChi2Cut);
}

//********************************************************************
void CoherentFitAlgorithm::ProcessMC(){
//********************************************************************

  std::cout << "--------------------" << std::endl;
  std::cout << "Processing MC coherent fit" << std::endl;
  std::cout << "--------------------" << std::endl;

  fmc->CreateCoherentSamples(fChi2Cut);
  fmc->GenerateTrueMCHistograms(fRMIN,fRMAX,fSTEP,fChi2Cut);
  fmc->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  fmc->SequentialCoherentFit();
}

//********************************************************************
void CoherentFitAlgorithm::ProcessData(){
//********************************************************************

  std::cout << "--------------------" << std::endl;
  std::cout << "Processing Data coherent fit" << std::endl;
  std::cout << "--------------------" << std::endl;

  fdata->CreateCoherentSamples(fChi2Cut);
  fdata->GenerateHistograms(fRMIN,fRMAX,fSTEP,fChi2Cut);
  fdata->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  fdata->DataCoherentFit(fmc);
}

