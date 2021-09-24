#include "CreateMiniTreeCosmics.hxx"
#include "Parameters.hxx"

//********************************************************************
CreateMiniTreeCosmics::CreateMiniTreeCosmics(int argc, char *argv[]):CreateMiniTreePD(argc, argv){
//********************************************************************

}

//********************************************************************
bool CreateMiniTreeCosmics::Initialize(){
//********************************************************************

  // call the base class
  CreateMiniTreePD::Initialize();
  
  _cutLength = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Length");
  _cutZmin   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Zmin");
  _cutZmax   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Zmax");
  
  return true;
}


//********************************************************************
bool CreateMiniTreeCosmics::CheckSaveParticle(const AnaParticleB& partB){
//********************************************************************

  const AnaParticlePD& part = *static_cast<const AnaParticlePD*>(&partB);

  if (part.ParentID==-1 &&
      part.Length > _cutLength &&
      part.PositionStart[2]>_cutZmin && part.PositionStart[2]<_cutZmax){

    return true;
  }
  else
    return false;
}

