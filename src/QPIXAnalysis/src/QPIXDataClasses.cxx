#define QPIXDataClasses_C

#include "QPIXDataClasses.hxx"
#include "AnalysisUtils.hxx"

//********************************************************************
AnaEventInfoQPIX::AnaEventInfoQPIX():AnaEventInfoPD(){
//********************************************************************

  IsBackground = false;
  HasPotasium  = false;
  HasGamma     = false;
  for(int i = 0; i < 7; i++){
    DetectedPhotons[i] = -999;
    PDSwf[i].clear();
  }
  SolarChain = unknown_chain;
  NuReaction = unknown_reaction;
}

//********************************************************************
AnaEventInfoQPIX::~AnaEventInfoQPIX(){
//********************************************************************

}

//********************************************************************
AnaEventInfoQPIX::AnaEventInfoQPIX(const AnaEventInfoQPIX& EventInfo):AnaEventInfoPD(EventInfo){
//********************************************************************

  IsBackground = EventInfo.IsBackground;
  HasPotasium  = EventInfo.HasPotasium;
  HasGamma     = EventInfo.HasGamma;
  for(int i = 0; i < 7; i++){
    DetectedPhotons[i] = EventInfo.DetectedPhotons[i];
    PDSwf[i]           = EventInfo.PDSwf[i];
  } 
  SolarChain   = EventInfo.SolarChain;
  NuReaction   = EventInfo.NuReaction;
}

//********************************************************************
void AnaEventInfoQPIX::Print() const{
//********************************************************************

  std::cout << "-------- AnaEventInfoQPIX --------- " << std::endl;
  AnaEventInfoPD::Print();
}

//********************************************************************
std::string AnaEventInfoQPIX::ChainToString(const int index) const{
//********************************************************************
  
  std::string result = "";
  if     (SolarChainEnum(index) == unknown_chain) result = "unknown_chain";
  else if(SolarChainEnum(index) == Be7)           result = "Be7";
  else if(SolarChainEnum(index) == pep)           result = "pep";
  else if(SolarChainEnum(index) == N13)           result = "N13";
  else if(SolarChainEnum(index) == O15)           result = "O15";
  else if(SolarChainEnum(index) == F17)           result = "F17";
  else if(SolarChainEnum(index) == B8)            result = "B8";
  else if(SolarChainEnum(index) == hep)           result = "hep";
  else                                            result = "unknown_chain";

  return result;
}

//********************************************************************
std::string AnaEventInfoQPIX::NuReactionToString(const int index) const{
//********************************************************************
  
  std::string result = "";
  if     (SolarNuReactionEnum(index) == unknown_reaction) result = "unknown_reaction";
  else if(SolarNuReactionEnum(index) == ES)               result = "ES";
  else if(SolarNuReactionEnum(index) == CC)               result = "CC";
  else                                                    result = "unknown_reaction";

  return result;
}
