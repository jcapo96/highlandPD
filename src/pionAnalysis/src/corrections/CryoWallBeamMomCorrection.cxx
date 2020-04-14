#include "CryoWallBeamMomCorrection.hxx"
#include "PionAnaDataClasses.hxx"
#include <cassert>

//#define DEBUG

//********************************************************************
CryoWallBeamMomCorrection::CryoWallBeamMomCorrection(){
//********************************************************************

}

//********************************************************************
void CryoWallBeamMomCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // Get beam
  AnaBeamPionAna* beam = static_cast<AnaBeamPionAna*>(spill.Beam);
  if (!beam)                        return; //?
    
  double corr = 0.0215;
      
  beam->BeamMomentumInTPC = beam->BeamMomentum - corr;
  
  //std::cout << "old momentum " << beam->BeamMomentum << " || new momentum " << beam->BeamMomentumInTPC << std::endl;

}


