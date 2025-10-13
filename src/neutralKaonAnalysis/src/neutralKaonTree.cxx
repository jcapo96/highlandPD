#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdDataClasses.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include "AnalysisUtils.hxx"
#include "ParticleId.hxx"
#include "HEPConstants.hxx"
#include <cmath>
#include <set>
#include <iostream>

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0(OutputManager& output, UInt_t nmax){
//********************************************************************

  // AddVarI(output, nk0, "number of neutral particle candidates");
  // AddVarI(output, k0nbrother, "number of K0 brothers");
  AddVarMaxSizeVI(output, k0hastrueobject, "K0 has true object", nk0, nmax);
  AddVarMaxSizeVI(output, k0hasrecoparticle, "K0 has reconstructed particle", nk0, nmax);
  AddVarMaxSizeVI(output, k0hasequivalenttrueobject, "K0 has equivalent true object", nk0, nmax);
  AddVarMaxSizeVI(output, k0id, "unique ID of neutral particle candidates", nk0, nmax);
  AddVarMaxSize3MF(output, k0recostartpos, "K0 reconstructed start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0truestartpos, "K0 true start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0recostartdir, "K0 reconstructed start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0truestartdir, "K0 true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0recoendpos, "K0 reconstructed end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0trueendpos, "K0 true end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0recoenddir, "K0 reconstructed end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0trueenddir, "K0 true end direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0truestartenddir, "K0 true start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0recostartenddir, "K0 reconstructed start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0fitstartenddir, "K0 fitted start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0recolength, "K0 reconstructed length", nk0, nmax);
  AddVarMaxSizeVF(output, k0truelength, "K0 true length", nk0, nmax);
  AddVarMaxSizeVF(output, k0truestartmom, "K0 true start momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0trueendmom, "K0 true end momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0recomass, "K0 reconstructed mass", nk0, nmax);
  AddVarMaxSizeVF(output, k0truemass, "K0 true mass", nk0, nmax);
  AddVarMaxSizeVI(output, k0truendau, "K0 true number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0trueproc, "K0 true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0trueendproc, "K0 true end process", nk0, nmax);
  AddVarMaxSizeVF(output, k0truerecodist, "K0 true-reco distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0impactparameter, "K0 impact parameter", nk0, nmax);
  AddVarMaxSizeVI(output, k0truepdg, "K0 true PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0recopdg, "K0 reconstructed PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0truegeneration, "K0 true generation", nk0, nmax);
  AddVarMaxSizeVI(output, k0nrecohits, "K0 number of reconstructed hits", nk0, nmax);

}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Par(OutputManager& output, UInt_t nmax){
  //********************************************************************

  // K0 parent variables
  AddVarMaxSize3MF(output, k0parrecostartpos, "K0 parent reconstructed start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0partruestartpos, "K0 parent true start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecoendpos, "K0 parent reconstructed end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0partrueendpos, "K0 parent true end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecostartdir, "K0 parent reconstructed start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0fitpardir, "K0 fitted parent direction from extrapolated track", nk0, nmax);
  AddVarMaxSize3MF(output, k0partruestartdir, "K0 parent true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecoenddir, "K0 parent reconstructed end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0partrueenddir, "K0 parent true end direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0parrecostartenddir, "K0 parent reconstructed start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruestartenddir, "K0 parent true start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0parrecolength, "K0 parent reconstructed length", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruelength, "K0 parent true length", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruestartmom, "K0 parent true start momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0partrueendmom, "K0 parent true end momentum", nk0, nmax);
  AddVarMaxSizeVI(output, k0partruendau, "K0 parent true number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0parrecondau, "K0 parent reconstructed number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0partruepdg, "K0 parent true PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0parrecopdg, "K0 parent reconstructed PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0partruepdgdau, "K0 parent true daughters PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0parrecopdgdau, "K0 parent reconstructed daughters PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0parisbeam, "K0 parent is beam particle", nk0, nmax);
  AddVarMaxSizeVI(output, k0beaminstpdg, "Beam instrumentation PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0partrueproc, "K0 parent true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0partrueendproc, "K0 parent true end process", nk0, nmax);
  AddVarMaxSizeVI(output, k0partruegeneration, "K0 parent true generation", nk0, nmax);

}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0vtxDaughter1(OutputManager& output, UInt_t nmax){
//********************************************************************

    // K0 daughter variables - COMMENTED OUT FOR NOW
    AddVarMaxSize3MF(output, k0dau1recostartpos, "K0 daughter1 reconstructed start position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1truestartpos, "K0 daughter1 true start position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1recoendpos, "K0 daughter1 reconstructed end position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1trueendpos, "K0 daughter1 true end position", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1recomom, "K0 daughter1 reconstructed momentum magnitude", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1recostartdir, "K0 daughter1 reconstructed start direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1fitdir, "K0 daughter1 fitted direction from vertex fit", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1truestartdir, "K0 daughter1 true start direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1recoenddir, "K0 daughter1 reconstructed end direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau1trueenddir, "K0 daughter1 true end direction", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1recolength, "K0 daughter1 reconstructed length", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1truelength, "K0 daughter1 true length", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1truestartmom, "K0 daughter1 true start momentum", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1trueendmom, "K0 daughter1 true end momentum", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1recondau, "K0 daughter1 reconstructed number of daughters", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1truendau, "K0 daughter1 true number of daughters", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1trueproc, "K0 daughter1 true process", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1trueendproc, "K0 daughter1 true end process", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1truepdg, "K0 daughter1 true PDG", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1parid, "K0 daughter1 parent ID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1istrack, "K0 daughter1 is track", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1chi2pion, "K0 daughter1 chi2 pion PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1nptchi2pion, "K0 daughter1 npt chi2 pion", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1chi2proton, "K0 daughter1 chi2 proton PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1nptchi2proton, "K0 daughter1 npt chi2 proton", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1chi2kaon, "K0 daughter1 chi2 kaon PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1nptchi2kaon, "K0 daughter1 npt chi2 kaon", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1avgdedx5cm, "K0 daughter1 avg dEdx RR<5cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1avgdedx10cm, "K0 daughter1 avg dEdx RR<10cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1avgdedx15cm, "K0 daughter1 avg dEdx RR<15cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1avgdedx25cm, "K0 daughter1 avg dEdx RR<25cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau1avgdedx50cm, "K0 daughter1 avg dEdx RR<50cm", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau1nhits, "K0 daughter1 number of hits in collection plane", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0vtxDaughter2(OutputManager& output, UInt_t nmax){
//********************************************************************

    AddVarMaxSize3MF(output, k0dau2recostartpos, "K0 daughter2 reconstructed start position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2truestartpos, "K0 daughter2 true start position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2recoendpos, "K0 daughter2 reconstructed end position", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2trueendpos, "K0 daughter2 true end position", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2recomom, "K0 daughter2 reconstructed momentum magnitude", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2recostartdir, "K0 daughter2 reconstructed start direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2fitdir, "K0 daughter2 fitted direction from vertex fit", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2truestartdir, "K0 daughter2 true start direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2recoenddir, "K0 daughter2 reconstructed end direction", nk0, nmax);
    AddVarMaxSize3MF(output, k0dau2trueenddir, "K0 daughter2 true end direction", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2recolength, "K0 daughter2 reconstructed length", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2truelength, "K0 daughter2 true length", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2truestartmom, "K0 daughter2 true start momentum", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2trueendmom, "K0 daughter2 true end momentum", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2recondau, "K0 daughter2 reconstructed number of daughters", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2truendau, "K0 daughter2 true number of daughters", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2trueproc, "K0 daughter2 true process", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2trueendproc, "K0 daughter2 true end process", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2truepdg, "K0 daughter2 true PDG", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2parid, "K0 daughter2 parent ID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2istrack, "K0 daughter2 is track", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2chi2pion, "K0 daughter2 chi2 pion PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2nptchi2pion, "K0 daughter2 npt chi2 pion", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2chi2proton, "K0 daughter2 chi2 proton PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2nptchi2proton, "K0 daughter2 npt chi2 proton", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2chi2kaon, "K0 daughter2 chi2 kaon PID", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2nptchi2kaon, "K0 daughter2 npt chi2 kaon", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2avgdedx5cm, "K0 daughter2 avg dEdx RR<5cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2avgdedx10cm, "K0 daughter2 avg dEdx RR<10cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2avgdedx15cm, "K0 daughter2 avg dEdx RR<15cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2avgdedx25cm, "K0 daughter2 avg dEdx RR<25cm", nk0, nmax);
    AddVarMaxSizeVF(output, k0dau2avgdedx50cm, "K0 daughter2 avg dEdx RR<50cm", nk0, nmax);
    AddVarMaxSizeVI(output, k0dau2nhits, "K0 daughter2 number of hits in collection plane", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Brother(OutputManager& output, UInt_t nmax){
//********************************************************************

   // Brother variables - using 2D matrix variables
   AddVarMaxSizeVI(output, k0reconbrother, "Number of K0 brothers", nk0, nmax);
   AddVarMaxSizeVI(output, k0truenbrother, "Number of K0 brothers", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconprot, "Number of K0 brothers that are reconstructed protons", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconpiplus, "Number of K0 brothers that are reconstructed pi+", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconpiminus, "Number of K0 brothers that are reconstructed pi-", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconprotchi2, "Number of K0 brothers compatible with proton (chi2/ndf<60)", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax){
//********************************************************************

  //Vertex system variables
  AddVarMaxSize3MF(output, k0vtxrecopos, "K0 vertex reconstructed position", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxtruepos, "K0 vertex true position", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxoriginaldistance, "K0 vertex original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueoriginaldistance, "K0 vertex true original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxminimumdistance, "K0 vertex minimum distance between fitted lines", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueminimumdistance, "K0 vertex true minimum distance between fitted lines", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxscore, "K0 vertex fit quality score", nk0, nmax);
  AddVarMaxSizeVF(output, k0neutralscore, "Neutral particle quality score (lower = more neutral-like)", nk0, nmax);
  AddVarMaxSizeVF(output, k0hitsalignment, "Alignment between hits in cylinder and neutral particle direction", nk0, nmax);
  AddVarMaxSizeVI(output, k0nhitsincylinder, "Number of hits in cylinder around neutral particle", nk0, nmax);
  AddVarMaxSizeVF(output, k0hitsavgdistance, "Average perpendicular distance of hits to neutral path (cm)", nk0, nmax);
  AddVarMaxSizeVF(output, k0hitsrmsdistance, "RMS of perpendicular distances of hits (cm)", nk0, nmax);
  AddVarMaxSizeVF(output, k0hitslongspan, "Longitudinal span fraction of hits along path", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpandorapos, "K0 vertex Pandora position", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdegbefore, "K0 vertex degeneracy before scoring", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdegafter, "K0 vertex degeneracy after scoring", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxnrecopart, "K0 vertex number of unique reco particles", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxrecomom, "Vertex system reconstructed momentum magnitude", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxtruemom, "Vertex system true momentum magnitude", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoenergy, "Vertex system reconstructed energy", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueenergy, "Vertex system true energy", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecomass, "Vertex system reconstructed mass", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtruemass, "Vertex system true mass", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxrecodir, "Vertex system reconstructed direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxtruedir, "Vertex system true direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoopening, "Vertex system reconstructed opening angle", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueopening, "Vertex system true opening angle", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoangle, "Vertex system reconstructed angle with beam", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueangle, "Vertex system true angle with beam", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxnpotpar, "Number of potential parent particles for vertex", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event, AnaBeamB* beam){
    // Fill all variables for all candidates
    if(candidate){
      neutralKaonTree::FillNeutralKaonVariables_K0(output, candidate);
      neutralKaonTree::FillNeutralKaonVariables_K0Par(output, candidate, beam);
      neutralKaonTree::FillNeutralKaonVariables_K0Brother(output, candidate->Parent);
      AnaVertexPD* vertex = candidate->Vertex;
      neutralKaonTree::FillNeutralKaonVariables_K0vtx(output, vertex);
      AnaParticlePD* daughter1Candidate = candidate->Vertex->Particles[0];
      neutralKaonTree::FillNeutralKaonVariables_K0vtxDaughter1(output, daughter1Candidate, vertex);
      AnaParticlePD* daughter2Candidate = candidate->Vertex->Particles[1];
      neutralKaonTree::FillNeutralKaonVariables_K0vtxDaughter2(output, daughter2Candidate, vertex);
    }
}
//********************************************************************

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************
    // Fill all variables for a single K0

  // K0 true variables (need to find associated true particle)
  Float_t k0truestartpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0truestartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueendpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0truestartenddir_val = -999.0;
  Float_t k0truelength_val = -999.0;
  Float_t k0truestartmom_val = -999.0;
  Float_t k0trueendmom_val = -999.0;
  Float_t k0truemass_val = -999.0;
  Int_t k0truendau_val = -999;
  Int_t k0trueproc_val = -999;
  Int_t k0trueendproc_val = -999;
  Int_t k0hastrueobject_val = -999;
  Int_t k0hasrecoparticle_val = -999;
  Int_t k0hasequivalenttrueobject_val = -999;
  Int_t k0truepdg_val = -999;
  Int_t k0recopdg_val = -999;
  Float_t k0truerecodist_val = -999.0;
  Int_t k0truegeneration_val = -999;

  if(candidate){
    // Fill unique ID of single candidate
    output.FillVectorVar(k0id, candidate->UniqueID);

    // K0 reconstructed start position (from AnaNeutralParticlePD inherited from AnaParticleB)
    Float_t k0recostartpos_val[3] = {-999.0, -999.0, -999.0};
    k0recostartpos_val[0] = candidate->PositionStart[0];
    k0recostartpos_val[1] = candidate->PositionStart[1];
    k0recostartpos_val[2] = candidate->PositionStart[2];
    output.FillMatrixVarFromArray(k0recostartpos, k0recostartpos_val, 3);

    //K0 reconstructed end position (vertex position)
    Float_t k0recoendpos_val[3] = {-999.0, -999.0, -999.0};
    if (candidate->Vertex) {
        k0recoendpos_val[0] = candidate->PositionEnd[0];
        k0recoendpos_val[1] = candidate->PositionEnd[1];
        k0recoendpos_val[2] = candidate->PositionEnd[2];
    }
    output.FillMatrixVarFromArray(k0recoendpos, k0recoendpos_val, 3);

    // K0 reconstructed start direction
    Float_t k0recostartdir_val[3] = {-999.0, -999.0, -999.0};
    k0recostartdir_val[0] = candidate->DirectionStart[0];
    k0recostartdir_val[1] = candidate->DirectionStart[1];
    k0recostartdir_val[2] = candidate->DirectionStart[2];
    output.FillMatrixVarFromArray(k0recostartdir, k0recostartdir_val, 3);

    // K0 reconstructed end direction
    Float_t k0recoenddir_val[3] = {-999.0, -999.0, -999.0};
    k0recoenddir_val[0] = candidate->DirectionEnd[0];
    k0recoenddir_val[1] = candidate->DirectionEnd[1];
    k0recoenddir_val[2] = candidate->DirectionEnd[2];
    output.FillMatrixVarFromArray(k0recoenddir, k0recoenddir_val, 3);

    // Calculate reconstructed scalar product
    Float_t k0recostartenddir_val = k0recostartdir_val[0]*k0recoenddir_val[0] + k0recostartdir_val[1]*k0recoenddir_val[1] + k0recostartdir_val[2]*k0recoenddir_val[2];
    output.FillVectorVar(k0recostartenddir, k0recostartenddir_val);

    // Calculate fitted start-end scalar product
    Float_t k0fitstartenddir_val = -999.0;
    if(candidate->FitParent && candidate->Vertex){
      k0fitstartenddir_val = candidate->FitParent->FitDirection[0]*candidate->Vertex->FitDirection[0] +
                             candidate->FitParent->FitDirection[1]*candidate->Vertex->FitDirection[1] +
                             candidate->FitParent->FitDirection[2]*candidate->Vertex->FitDirection[2];
    }
    output.FillVectorVar(k0fitstartenddir, k0fitstartenddir_val);

    // K0 reconstructed length
    Float_t k0recolength_val = -999.0;
    k0recolength_val = sqrt(pow(k0recoendpos_val[0] - k0recostartpos_val[0], 2) + pow(k0recoendpos_val[1] - k0recostartpos_val[1], 2) + pow(k0recoendpos_val[2] - k0recostartpos_val[2], 2));
    output.FillVectorVar(k0recolength, k0recolength_val);

    // K0 reconstructed mass
    Float_t k0recomass_val = -999.0;
    // Calculate invariant mass from vertex daughters assuming pion hypothesis
    if(candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
      AnaParticlePD* daughter1 = candidate->Vertex->Particles[0];
      AnaParticlePD* daughter2 = candidate->Vertex->Particles[1];

      if(daughter1 && daughter2) {
        // Use framework mass constant (in GeV)
        const double pionMass = 0.13957;

        // Get reconstructed momenta (in GeV/c from EnsureParticleMomentum)
        double p1 = daughter1->Momentum;
        double p2 = daughter2->Momentum;

        // Calculate energies assuming pion mass (all in GeV)
        double E1 = sqrt(p1*p1 + pionMass*pionMass);
        double E2 = sqrt(p2*p2 + pionMass*pionMass);

        // Get momentum vectors (in GeV/c)
        TVector3 p1vec(daughter1->DirectionStart[0] * p1,
                      daughter1->DirectionStart[1] * p1,
                      daughter1->DirectionStart[2] * p1);
        TVector3 p2vec(daughter2->DirectionStart[0] * p2,
                      daughter2->DirectionStart[1] * p2,
                      daughter2->DirectionStart[2] * p2);

        // Calculate invariant mass (in MeV/c^2)
        TVector3 pTot = p1vec + p2vec;
        double ETot = E1 + E2;
        double massGeV = sqrt(ETot*ETot - pTot.Mag2());

        // Convert to MeV/c^2
        k0recomass_val = massGeV * 1000;
      }
    }
    output.FillVectorVar(k0recomass, k0recomass_val);

    // K0 impact parameter
    Float_t k0impactparameter_val = candidate->ImpactParameter;
    output.FillVectorVar(k0impactparameter, k0impactparameter_val);

    // K0 number of reconstructed hits in cylinder
    Int_t k0nrecohits_val = candidate->NRecoHitsInVertex;
    output.FillVectorVar(k0nrecohits, k0nrecohits_val);

    // K0 neutral particle quality score
    Float_t k0neutralscore_val = candidate->NeutralScore;
    output.FillVectorVar(k0neutralscore, k0neutralscore_val);

    // K0 hits alignment (dot product between hits direction and neutral particle direction)
    Float_t k0hitsalignment_val = candidate->HitsAlignment;
    output.FillVectorVar(k0hitsalignment, k0hitsalignment_val);

    // K0 number of hits in cylinder
    Int_t k0nhitsincylinder_val = candidate->NHitsInCylinder;
    output.FillVectorVar(k0nhitsincylinder, k0nhitsincylinder_val);

    // K0 average perpendicular distance of hits
    Float_t k0hitsavgdistance_val = candidate->HitsAvgDistance;
    output.FillVectorVar(k0hitsavgdistance, k0hitsavgdistance_val);

    // K0 RMS of perpendicular distances
    Float_t k0hitsrmsdistance_val = candidate->HitsRMSDistance;
    output.FillVectorVar(k0hitsrmsdistance, k0hitsrmsdistance_val);

    // K0 longitudinal span fraction
    Float_t k0hitslongspan_val = candidate->HitsLongitudinalSpan;
    output.FillVectorVar(k0hitslongspan, k0hitslongspan_val);

    // Determine if we have an equivalent true object
    if(candidate->TrueEquivalentNeutralParticle){
      k0hasequivalenttrueobject_val = 1;
    }
    else{
      k0hasequivalenttrueobject_val = 0;
    }


    // Determine if we have a true object
    if(candidate->TrueObject){
        k0hastrueobject_val = 1;
        if(candidate->RecoParticle){
          k0hasrecoparticle_val = 1;
        }
        else{
          k0hasrecoparticle_val = 0;
        }
        // std::cout << "DEBUG: Found TrueObject" << std::endl;
        AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(candidate->TrueObject);
          if(trueNeutralParticle){
            k0truestartpos_val[0] = trueNeutralParticle->Position[0];
            k0truestartpos_val[1] = trueNeutralParticle->Position[1];
            k0truestartpos_val[2] = trueNeutralParticle->Position[2];
            k0truestartdir_val[0] = trueNeutralParticle->Direction[0];
            k0truestartdir_val[1] = trueNeutralParticle->Direction[1];
            k0truestartdir_val[2] = trueNeutralParticle->Direction[2];
            k0trueendpos_val[0] = trueNeutralParticle->PositionEnd[0];
            k0trueendpos_val[1] = trueNeutralParticle->PositionEnd[1];
            k0trueendpos_val[2] = trueNeutralParticle->PositionEnd[2];
            k0trueenddir_val[0] = trueNeutralParticle->DirectionEnd[0];
            k0trueenddir_val[1] = trueNeutralParticle->DirectionEnd[1];
            k0trueenddir_val[2] = trueNeutralParticle->DirectionEnd[2];
            k0truelength_val = sqrt(pow(k0trueendpos_val[0] - k0truestartpos_val[0], 2) + pow(k0trueendpos_val[1] - k0truestartpos_val[1], 2) + pow(k0trueendpos_val[2] - k0truestartpos_val[2], 2));
            k0truestartmom_val = trueNeutralParticle->Momentum;
          k0trueendmom_val = trueNeutralParticle->Momentum;
            k0truendau_val = trueNeutralParticle->Daughters.size();
            k0trueproc_val = static_cast<Int_t>(trueNeutralParticle->ProcessStart);
            k0trueendproc_val = static_cast<Int_t>(trueNeutralParticle->ProcessEnd);
            k0truestartenddir_val = trueNeutralParticle->Direction[0]*trueNeutralParticle->DirectionEnd[0] + trueNeutralParticle->Direction[1]*trueNeutralParticle->DirectionEnd[1] + trueNeutralParticle->Direction[2]*trueNeutralParticle->DirectionEnd[2];
            k0truepdg_val = trueNeutralParticle->PDG;
          k0truerecodist_val = sqrt(pow(k0trueendpos_val[0] - k0recoendpos_val[0], 2) + pow(k0trueendpos_val[1] - k0recoendpos_val[1], 2) + pow(k0trueendpos_val[2] - k0recoendpos_val[2], 2));
          k0truegeneration_val = trueNeutralParticle->Generation;

          // Calculate true mass from vertex daughters assuming pion hypothesis
          if(candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
            AnaParticlePD* daughter1 = candidate->Vertex->Particles[0];
            AnaParticlePD* daughter2 = candidate->Vertex->Particles[1];

            if(daughter1 && daughter2 && daughter1->TrueObject && daughter2->TrueObject) {
              AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
              AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

              if(trueDaughter1 && trueDaughter2) {
                // Use framework mass constant (in GeV)
                const double pionMass = 0.13957;

                // Get momenta (in framework base units = GeV/c)
                double p1 = trueDaughter1->Momentum;
                double p2 = trueDaughter2->Momentum;

                // std::cout << "DEBUG: True daughter 1 momentum: " << p1 << std::endl;
                // std::cout << "DEBUG: True daughter 2 momentum: " << p2 << std::endl;

                // Calculate energies assuming pion mass (all in MeV)
                double E1 = sqrt(p1*p1 + pionMass*pionMass);
                double E2 = sqrt(p2*p2 + pionMass*pionMass);

                // std::cout << "DEBUG: True daughter 1 energy: " << E1 << std::endl;
                // std::cout << "DEBUG: True daughter 2 energy: " << E2 << std::endl;

                // Get momentum vectors (in GeV/c)
                TVector3 p1vec(trueDaughter1->Direction[0] * p1,
                              trueDaughter1->Direction[1] * p1,
                              trueDaughter1->Direction[2] * p1);
                TVector3 p2vec(trueDaughter2->Direction[0] * p2,
                              trueDaughter2->Direction[1] * p2,
                              trueDaughter2->Direction[2] * p2);

                // Calculate invariant mass (in GeV/c^2)
                TVector3 pTot = p1vec + p2vec;
                double ETot = E1 + E2;
                // std::cout << "DEBUG: True daughter 1 total energy: " << ETot << std::endl;
                double massGeV = sqrt(ETot*ETot - pTot.Mag2());

                // std::cout << "DEBUG: True daughter 1 invariant mass: " << massGeV << std::endl;

                // Convert to GeV/c^2
                k0truemass_val = massGeV * 1000; // Convert to MeV/c^2
                // std::cout << "DEBUG: True daughter 1 invariant mass: " << k0truemass_val << std::endl;
              }
            }
          }

          // Get the reco object associated to the true object to see if we have a reconstructed PDG
          if(trueNeutralParticle->ReconParticles.size() > 0){
            AnaParticlePD* recoCandidate = static_cast<AnaParticlePD*>(trueNeutralParticle->ReconParticles[0]);
            if(recoCandidate){
              if(recoCandidate->ReconPDG[0] != -999){
                k0recopdg_val = recoCandidate->ReconPDG[0];
              }
              else if(recoCandidate->ReconPDG[1] != -999){
                k0recopdg_val = recoCandidate->ReconPDG[1];
              }
              else if(recoCandidate->ReconPDG[2] != -999){
                k0recopdg_val = recoCandidate->ReconPDG[2];
            }
            else{
                k0recopdg_val = -999;
              }
            }
          }
        }
      }
      // If no true object, check if there is an equivalent true object
      else{
        // k0hastrueobject_val remains 0 (has equivalent but not direct true object)
        k0hastrueobject_val = 0;
        AnaTrueEquivalentNeutralParticlePD* trueEquivalentNeutralParticle = static_cast<AnaTrueEquivalentNeutralParticlePD*>(candidate->TrueEquivalentNeutralParticle);
        if(trueEquivalentNeutralParticle){
          k0truestartpos_val[0] = trueEquivalentNeutralParticle->Position[0];
          k0truestartpos_val[1] = trueEquivalentNeutralParticle->Position[1];
          k0truestartpos_val[2] = trueEquivalentNeutralParticle->Position[2];
          k0truestartdir_val[0] = trueEquivalentNeutralParticle->Direction[0];
          k0truestartdir_val[1] = trueEquivalentNeutralParticle->Direction[1];
          k0truestartdir_val[2] = trueEquivalentNeutralParticle->Direction[2];
          k0trueendpos_val[0] = trueEquivalentNeutralParticle->PositionEnd[0];
          k0trueendpos_val[1] = trueEquivalentNeutralParticle->PositionEnd[1];
          k0trueendpos_val[2] = trueEquivalentNeutralParticle->PositionEnd[2];
          k0trueenddir_val[0] = trueEquivalentNeutralParticle->DirectionEnd[0];
          k0trueenddir_val[1] = trueEquivalentNeutralParticle->DirectionEnd[1];
          k0trueenddir_val[2] = trueEquivalentNeutralParticle->DirectionEnd[2];
          k0truelength_val = sqrt(pow(k0trueendpos_val[0] - k0truestartpos_val[0], 2) + pow(k0trueendpos_val[1] - k0truestartpos_val[1], 2) + pow(k0trueendpos_val[2] - k0truestartpos_val[2], 2));
          k0truestartmom_val = trueEquivalentNeutralParticle->Momentum;
          k0trueendmom_val = trueEquivalentNeutralParticle->Momentum;
          k0truendau_val = trueEquivalentNeutralParticle->TrueEquivalentVertex->TrueParticles.size();
          k0trueproc_val = trueEquivalentNeutralParticle->Process;
          // Note: AnaTrueEquivalentNeutralParticlePD only has Process, not ProcessStart/ProcessEnd
          // k0trueendproc_val remains -999 for equivalent true particles
          k0truestartenddir_val = trueEquivalentNeutralParticle->Direction[0]*trueEquivalentNeutralParticle->DirectionEnd[0] + trueEquivalentNeutralParticle->Direction[1]*trueEquivalentNeutralParticle->DirectionEnd[1] + trueEquivalentNeutralParticle->Direction[2]*trueEquivalentNeutralParticle->DirectionEnd[2];
          k0truerecodist_val = sqrt(pow(k0trueendpos_val[0] - k0recoendpos_val[0], 2) + pow(k0trueendpos_val[1] - k0recoendpos_val[1], 2) + pow(k0trueendpos_val[2] - k0recoendpos_val[2], 2));
          k0truepdg_val = trueEquivalentNeutralParticle->PDG;

          // Calculate true mass from equivalent vertex true particles assuming pion hypothesis
          if(trueEquivalentNeutralParticle->TrueEquivalentVertex &&
             trueEquivalentNeutralParticle->TrueEquivalentVertex->TrueParticles.size() >= 2) {
            AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(trueEquivalentNeutralParticle->TrueEquivalentVertex->TrueParticles[0]);
            AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(trueEquivalentNeutralParticle->TrueEquivalentVertex->TrueParticles[1]);

            if(trueDaughter1 && trueDaughter2) {
              // Use framework mass constant (in GeV)
              const double pionMass = 0.13957;

              // Get momenta (in framework base units = GeV/c)
              double p1 = trueDaughter1->Momentum;
              double p2 = trueDaughter2->Momentum;

              // Calculate energies assuming pion mass (all in MeV)
              double E1 = sqrt(p1*p1 + pionMass*pionMass);
              double E2 = sqrt(p2*p2 + pionMass*pionMass);

              // Get momentum vectors (in GeV/c)
              TVector3 p1vec(trueDaughter1->Direction[0] * p1,
                            trueDaughter1->Direction[1] * p1,
                            trueDaughter1->Direction[2] * p1);
              TVector3 p2vec(trueDaughter2->Direction[0] * p2,
                            trueDaughter2->Direction[1] * p2,
                            trueDaughter2->Direction[2] * p2);

              // Calculate invariant mass (in MeV/c^2)
              TVector3 pTot = p1vec + p2vec;
              double ETot = E1 + E2;
              double massGeV = sqrt(ETot*ETot - pTot.Mag2());

              // Convert to GeV/c^2
              k0truemass_val = massGeV * 1000; // Convert to MeV/c^2
            }
          }
        }
      }
    }

  // Handle null candidate case - fill with default values
      output.FillMatrixVarFromArray(k0truestartpos, k0truestartpos_val, 3);
      output.FillMatrixVarFromArray(k0truestartdir, k0truestartdir_val, 3);
      output.FillMatrixVarFromArray(k0trueendpos, k0trueendpos_val, 3);
      output.FillMatrixVarFromArray(k0trueenddir, k0trueenddir_val, 3);
      output.FillVectorVar(k0truestartenddir, k0truestartenddir_val);
      output.FillVectorVar(k0truelength, k0truelength_val);
      output.FillVectorVar(k0truestartmom, k0truestartmom_val);
      output.FillVectorVar(k0trueendmom, k0trueendmom_val);
      output.FillVectorVar(k0truemass, k0truemass_val);
      output.FillVectorVar(k0truendau, k0truendau_val);
      output.FillVectorVar(k0trueproc, k0trueproc_val);
      output.FillVectorVar(k0trueendproc, k0trueendproc_val);
      output.FillVectorVar(k0hastrueobject, k0hastrueobject_val);
      output.FillVectorVar(k0hasrecoparticle, k0hasrecoparticle_val);
      output.FillVectorVar(k0hasequivalenttrueobject, k0hasequivalenttrueobject_val);
      output.FillVectorVar(k0truepdg, k0truepdg_val);
      output.FillVectorVar(k0recopdg, k0recopdg_val);
      output.FillVectorVar(k0truerecodist, k0truerecodist_val);
      output.FillVectorVar(k0truegeneration, k0truegeneration_val);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Par(OutputManager& output, AnaNeutralParticlePD* neutralCandidate, AnaBeamB* beam){
//********************************************************************
    // Fill all variables for a single K0 parent
    AnaParticlePD* parentCandidate = neutralCandidate ? neutralCandidate->Parent : NULL;
    if(parentCandidate){

      Float_t k0parrecostartpos_val[3] = {-999.0, -999.0, -999.0};
      k0parrecostartpos_val[0] = parentCandidate->PositionStart[0];
      k0parrecostartpos_val[1] = parentCandidate->PositionStart[1];
      k0parrecostartpos_val[2] = parentCandidate->PositionStart[2];
      output.FillMatrixVarFromArray(k0parrecostartpos, k0parrecostartpos_val, 3);

      Float_t k0parrecoendpos_val[3] = {-999.0, -999.0, -999.0};
      k0parrecoendpos_val[0] = parentCandidate->PositionEnd[0];
      k0parrecoendpos_val[1] = parentCandidate->PositionEnd[1];
      k0parrecoendpos_val[2] = parentCandidate->PositionEnd[2];
      output.FillMatrixVarFromArray(k0parrecoendpos, k0parrecoendpos_val, 3);

      Float_t k0parrecostartdir_val[3] = {-999.0, -999.0, -999.0};
      k0parrecostartdir_val[0] = parentCandidate->DirectionStart[0];
      k0parrecostartdir_val[1] = parentCandidate->DirectionStart[1];
      k0parrecostartdir_val[2] = parentCandidate->DirectionStart[2];
      output.FillMatrixVarFromArray(k0parrecostartdir, k0parrecostartdir_val, 3);

      // K0 fitted parent direction from extrapolated track
      Float_t k0fitpardir_val[3] = {-999.0, -999.0, -999.0};
      if(neutralCandidate && neutralCandidate->FitParent){
        k0fitpardir_val[0] = neutralCandidate->FitParent->FitDirection[0];
        k0fitpardir_val[1] = neutralCandidate->FitParent->FitDirection[1];
        k0fitpardir_val[2] = neutralCandidate->FitParent->FitDirection[2];
      }
      output.FillMatrixVarFromArray(k0fitpardir, k0fitpardir_val, 3);

      Float_t k0parrecoenddir_val[3] = {-999.0, -999.0, -999.0};
      k0parrecoenddir_val[0] = parentCandidate->DirectionEnd[0];
      k0parrecoenddir_val[1] = parentCandidate->DirectionEnd[1];
      k0parrecoenddir_val[2] = parentCandidate->DirectionEnd[2];
      output.FillMatrixVarFromArray(k0parrecoenddir, k0parrecoenddir_val, 3);

      // Calculate reconstructed start-end direction scalar product
      Float_t k0parrecostartenddir_val = k0parrecostartdir_val[0]*k0parrecoenddir_val[0] + k0parrecostartdir_val[1]*k0parrecoenddir_val[1] + k0parrecostartdir_val[2]*k0parrecoenddir_val[2];
      output.FillVectorVar(k0parrecostartenddir, k0parrecostartenddir_val);

      Float_t k0parrecolength_val = -999.0;
      k0parrecolength_val = parentCandidate->Length;
      output.FillVectorVar(k0parrecolength, k0parrecolength_val);

      Float_t k0partruestartmom_val = -999.0;
      k0partruestartmom_val = parentCandidate->Momentum;
      output.FillVectorVar(k0partruestartmom, k0partruestartmom_val);

      Float_t k0partrueendmom_val = -999.0;
      k0partrueendmom_val = parentCandidate->MomentumEnd;
      output.FillVectorVar(k0partrueendmom, k0partrueendmom_val);

    Int_t k0parrecondau_val = -999.0;
      k0parrecondau_val = parentCandidate->Daughters.size();
      output.FillVectorVar(k0parrecondau, k0parrecondau_val);

    Int_t k0parrecopdg_val = -999;
    if(parentCandidate->ReconPDG[0] != -999){
      k0parrecopdg_val = parentCandidate->ReconPDG[0];
    }
    else if(parentCandidate->ReconPDG[1] != -999){
      k0parrecopdg_val = parentCandidate->ReconPDG[1];
    }
    else if(parentCandidate->ReconPDG[2] != -999){
      k0parrecopdg_val = parentCandidate->ReconPDG[2];
    }
    else{
      k0parrecopdg_val = -999;
    }
      output.FillVectorVar(k0parrecopdg, k0parrecopdg_val);

    Int_t k0parisbeam_val = -999;
    k0parisbeam_val = parentCandidate->isPandora;
      output.FillVectorVar(k0parisbeam, k0parisbeam_val);

    // Fill beam instrumentation PDG
    Int_t k0beaminstpdg_val = -999;
    if(beam){
      AnaBeamPD* beamPD = static_cast<AnaBeamPD*>(beam);
      if(beamPD && (int)beamPD->PDGs.size() > 0){
        k0beaminstpdg_val = beamPD->PDGs[0];
      }
    }
    output.FillVectorVar(k0beaminstpdg, k0beaminstpdg_val);

    Int_t k0partruegeneration_val = -999;

      Float_t k0partruestartpos_val[3] = {-999.0, -999.0, -999.0};
      Float_t k0partruestartdir_val[3] = {-999.0, -999.0, -999.0};
      Float_t k0partrueendpos_val[3] = {-999.0, -999.0, -999.0};
      Float_t k0partrueenddir_val[3] = {-999.0, -999.0, -999.0};
      Float_t k0partruestartenddir_val = -999.0;
      Float_t k0partruelength_val = -999.0;
      Int_t k0partruendau_val = -999;
      Int_t k0partrueproc_val = -999;
      Int_t k0partrueendproc_val = -999;
      Int_t k0partruepdg_val = -999;

      AnaTrueParticlePD* trueParentCandidate = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);
      if(trueParentCandidate){
        k0partruestartpos_val[0] = trueParentCandidate->Position[0];
        k0partruestartpos_val[1] = trueParentCandidate->Position[1];
        k0partruestartpos_val[2] = trueParentCandidate->Position[2];
        k0partruestartdir_val[0] = trueParentCandidate->Direction[0];
        k0partruestartdir_val[1] = trueParentCandidate->Direction[1];
        k0partruestartdir_val[2] = trueParentCandidate->Direction[2];
      k0partrueendpos_val[0] = trueParentCandidate->PositionEnd[0];
      k0partrueendpos_val[1] = trueParentCandidate->PositionEnd[1];
      k0partrueendpos_val[2] = trueParentCandidate->PositionEnd[2];
      k0partrueenddir_val[0] = trueParentCandidate->DirectionEnd[0];
      k0partrueenddir_val[1] = trueParentCandidate->DirectionEnd[1];
      k0partrueenddir_val[2] = trueParentCandidate->DirectionEnd[2];
      k0partruelength_val = sqrt(pow(k0partrueendpos_val[0] - k0partruestartpos_val[0], 2) + pow(k0partrueendpos_val[1] - k0partruestartpos_val[1], 2) + pow(k0partrueendpos_val[2] - k0partruestartpos_val[2], 2));
      k0partruestartenddir_val = k0partruestartdir_val[0]*k0partrueenddir_val[0] + k0partruestartdir_val[1]*k0partrueenddir_val[1] + k0partruestartdir_val[2]*k0partrueenddir_val[2];
      k0partruestartmom_val = trueParentCandidate->Momentum;
      k0partrueendmom_val = trueParentCandidate->Momentum;
      k0partruendau_val = trueParentCandidate->Daughters.size();
      k0partrueproc_val = static_cast<Int_t>(trueParentCandidate->ProcessStart);
      k0partrueendproc_val = static_cast<Int_t>(trueParentCandidate->ProcessEnd);
      k0partruepdg_val = trueParentCandidate->PDG;
      k0partruegeneration_val = trueParentCandidate->Generation;
  }

  output.FillMatrixVarFromArray(k0partruestartpos, k0partruestartpos_val, 3);
  output.FillMatrixVarFromArray(k0partruestartdir, k0partruestartdir_val, 3);
  output.FillMatrixVarFromArray(k0partrueendpos, k0partrueendpos_val, 3);
  output.FillMatrixVarFromArray(k0partrueenddir, k0partrueenddir_val, 3);
    output.FillVectorVar(k0partruestartenddir, k0partruestartenddir_val);
  output.FillVectorVar(k0partruelength, k0partruelength_val);
  output.FillVectorVar(k0partruestartmom, k0partruestartmom_val);
  output.FillVectorVar(k0partrueendmom, k0partrueendmom_val);
  output.FillVectorVar(k0partruendau, k0partruendau_val);
  output.FillVectorVar(k0partrueproc, k0partrueproc_val);
  output.FillVectorVar(k0partrueendproc, k0partrueendproc_val);
  output.FillVectorVar(k0partruepdg, k0partruepdg_val);
  output.FillVectorVar(k0partruegeneration, k0partruegeneration_val);
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaVertexPD* vertex){
    // Fill all variables for a single K0 vertex
  if(vertex){
    Float_t k0vtxrecopos_val[3] = {-999.0, -999.0, -999.0};
    k0vtxrecopos_val[0] = vertex->Position[0];
    k0vtxrecopos_val[1] = vertex->Position[1];
    k0vtxrecopos_val[2] = vertex->Position[2];

    // Get vertex momentum, mass, energy, and angle from vertex attributes
    // These are now calculated when the vertex is created in pdNeutralUtils
    Float_t k0vtxrecomom_val[3] = {vertex->Momentum[0], vertex->Momentum[1], vertex->Momentum[2]};
    Float_t k0vtxrecoenergy_val = vertex->Energy;

    // Calculate reconstructed mass from vertex particles
    Float_t k0vtxrecomass_val = -999.0;
    if(vertex->Particles.size() >= 2) {
      AnaParticlePD* particle1 = vertex->Particles[0];
      AnaParticlePD* particle2 = vertex->Particles[1];

      if(particle1 && particle2 &&
         particle1->Momentum > 0 && particle2->Momentum > 0 &&
         particle1->Momentum != -999 && particle2->Momentum != -999) {
        // Calculate invariant mass
        const double pionMass = 0.13957; // GeV
        TVector3 p1vec(particle1->DirectionStart[0] * particle1->Momentum,
                      particle1->DirectionStart[1] * particle1->Momentum,
                      particle1->DirectionStart[2] * particle1->Momentum);
        TVector3 p2vec(particle2->DirectionStart[0] * particle2->Momentum,
                      particle2->DirectionStart[1] * particle2->Momentum,
                      particle2->DirectionStart[2] * particle2->Momentum);
        TVector3 pTot = p1vec + p2vec;
        double E1 = sqrt(particle1->Momentum*particle1->Momentum + pionMass*pionMass);
        double E2 = sqrt(particle2->Momentum*particle2->Momentum + pionMass*pionMass);
        double ETot = E1 + E2;
        double massGeV = sqrt(ETot*ETot - pTot.Mag2());
        k0vtxrecomass_val = massGeV * 1000; // Convert to MeV/c^2
      }
    }

    Float_t k0vtxrecodir_val[3] = {vertex->Direction[0], vertex->Direction[1], vertex->Direction[2]};
    Float_t k0vtxrecoopening_val = vertex->OpeningAngle;
    Float_t k0vtxrecoangle_val = vertex->AngleWithBeam;

    Float_t k0vtxoriginaldistance_val = -999.0;
    k0vtxoriginaldistance_val = vertex->OriginalDistance;

    Float_t k0vtxminimumdistance_val = -999.0;
    k0vtxminimumdistance_val = vertex->MinimumDistance;

    Float_t k0vtxscore_val = -999.0;
    k0vtxscore_val = vertex->Score;

    Int_t k0vtxnpotpar_val = -999;
    k0vtxnpotpar_val = vertex->NPotentialParents;

    Float_t k0vtxtruepos_val[3] = {-999.0, -999.0, -999.0};
    Float_t k0vtxtruedir_val[3] = {-999.0, -999.0, -999.0};
    Float_t k0vtxtruemom_val[3] = {-999.0, -999.0, -999.0};
    Float_t k0vtxtrueenergy_val = -999.0;
    Float_t k0vtxtruemass_val = -999.0;
    Float_t k0vtxtrueoriginaldistance_val = -999.0;
    Float_t k0vtxtrueminimumdistance_val = -999.0;
    Float_t k0vtxtrueopening_val = -999.0;
    Float_t k0vtxtrueangle_val = -999.0;

    if(vertex->TrueEquivalentVertex){
      AnaTrueEquivalentVertexPD* trueEquivalentVertex = static_cast<AnaTrueEquivalentVertexPD*>(vertex->TrueEquivalentVertex);
      if(trueEquivalentVertex){
        k0vtxtruepos_val[0] = trueEquivalentVertex->Position[0];
        k0vtxtruepos_val[1] = trueEquivalentVertex->Position[1];
        k0vtxtruepos_val[2] = trueEquivalentVertex->Position[2];
        k0vtxtruedir_val[0] = trueEquivalentVertex->Direction[0];
        k0vtxtruedir_val[1] = trueEquivalentVertex->Direction[1];
        k0vtxtruedir_val[2] = trueEquivalentVertex->Direction[2];
        k0vtxtrueoriginaldistance_val = trueEquivalentVertex->OriginalDistance;
        k0vtxtrueminimumdistance_val = trueEquivalentVertex->MinimumDistance;
        k0vtxtrueopening_val = trueEquivalentVertex->OpeningAngle;

        // Calculate true angle with beam
        if(k0vtxtruedir_val[0] != -999.0) {
          TVector3 vtxTrueDir(k0vtxtruedir_val[0], k0vtxtruedir_val[1], k0vtxtruedir_val[2]);
          TVector3 beamDir(0, 0, 1); // Beam along +z
          k0vtxtrueangle_val = vtxTrueDir.Angle(beamDir);
        }

        // Calculate true momentum, energy, and mass if we have two true particles
        if(trueEquivalentVertex->TrueParticles.size() >= 2) {
          AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(trueEquivalentVertex->TrueParticles[0]);
          AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(trueEquivalentVertex->TrueParticles[1]);

          if(trueParticle1 && trueParticle2) {
            const double pionMass = 0.13957; // GeV

            // Calculate total momentum vector
            TVector3 p1vec(trueParticle1->Direction[0] * trueParticle1->Momentum,
                          trueParticle1->Direction[1] * trueParticle1->Momentum,
                          trueParticle1->Direction[2] * trueParticle1->Momentum);
            TVector3 p2vec(trueParticle2->Direction[0] * trueParticle2->Momentum,
                          trueParticle2->Direction[1] * trueParticle2->Momentum,
                          trueParticle2->Direction[2] * trueParticle2->Momentum);
            TVector3 pTot = p1vec + p2vec;

            k0vtxtruemom_val[0] = pTot.X();
            k0vtxtruemom_val[1] = pTot.Y();
            k0vtxtruemom_val[2] = pTot.Z();

            // Calculate total energy
            double E1 = sqrt(trueParticle1->Momentum*trueParticle1->Momentum + pionMass*pionMass);
            double E2 = sqrt(trueParticle2->Momentum*trueParticle2->Momentum + pionMass*pionMass);
            double ETot = E1 + E2;
            k0vtxtrueenergy_val = ETot * 1000; // Convert to MeV

            // Calculate invariant mass
            double massGeV = sqrt(ETot*ETot - pTot.Mag2());
            k0vtxtruemass_val = massGeV * 1000; // Convert to MeV/c^2
          }
        }
      }
    }

    output.FillMatrixVarFromArray(k0vtxrecopos, k0vtxrecopos_val, 3);
    output.FillMatrixVarFromArray(k0vtxrecomom, k0vtxrecomom_val, 3);
    output.FillVectorVar(k0vtxrecoenergy, k0vtxrecoenergy_val);
    output.FillMatrixVarFromArray(k0vtxrecodir, k0vtxrecodir_val, 3);
    output.FillVectorVar(k0vtxrecomass, k0vtxrecomass_val);
    output.FillVectorVar(k0vtxrecoopening, k0vtxrecoopening_val);
    output.FillVectorVar(k0vtxrecoangle, k0vtxrecoangle_val);
    output.FillVectorVar(k0vtxoriginaldistance, k0vtxoriginaldistance_val);
    output.FillVectorVar(k0vtxminimumdistance, k0vtxminimumdistance_val);
    output.FillVectorVar(k0vtxscore, k0vtxscore_val);
    output.FillVectorVar(k0vtxnpotpar, k0vtxnpotpar_val);

    output.FillMatrixVarFromArray(k0vtxtruepos, k0vtxtruepos_val, 3);
    output.FillMatrixVarFromArray(k0vtxtruemom, k0vtxtruemom_val, 3);
    output.FillVectorVar(k0vtxtrueenergy, k0vtxtrueenergy_val);
    output.FillMatrixVarFromArray(k0vtxtruedir, k0vtxtruedir_val, 3);
    output.FillVectorVar(k0vtxtruemass, k0vtxtruemass_val);
    output.FillVectorVar(k0vtxtrueoriginaldistance, k0vtxtrueoriginaldistance_val);
    output.FillVectorVar(k0vtxtrueminimumdistance, k0vtxtrueminimumdistance_val);
    output.FillVectorVar(k0vtxtrueopening, k0vtxtrueopening_val);
    output.FillVectorVar(k0vtxtrueangle, k0vtxtrueangle_val);

    // Fill Pandora position
    output.FillMatrixVarFromArray(k0vtxpandorapos, vertex->PositionPandora, 3);

    // Fill degeneracy variables
    output.FillVectorVar(k0vtxdegbefore, vertex->DegeneracyBeforeScoring);
    output.FillVectorVar(k0vtxdegafter, vertex->DegeneracyAfterScoring);
    output.FillVectorVar(k0vtxnrecopart, vertex->NRecoParticles);
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Brother(OutputManager& output, AnaParticlePD* parentCandidate){
    // Fill all variables for a single K0 brother
  if(parentCandidate){
    output.FillVectorVar(k0reconbrother, (Int_t)parentCandidate->Daughters.size());
    AnaTrueParticlePD* trueParentCandidate = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);

    // Get minimum length parameter for counting brothers
    double minBrotherLength = ND::params().GetParameterD("neutralKaonAnalysis.MinBrotherLength");

    UInt_t nReconProtons = 0;
    UInt_t nRecoPiPlus = 0;
    UInt_t nRecoPiMinus = 0;
    UInt_t nReconProtChi2 = 0;

    for(size_t i = 0; i < parentCandidate->Daughters.size(); i++){
      AnaParticlePD* daughterCandidate = static_cast<AnaParticlePD*>(parentCandidate->Daughters[i]);
      if(daughterCandidate){
        // Check if daughter is long enough
        if(daughterCandidate->Length < minBrotherLength){
          continue; // Skip short daughters
        }

        // Count by true PDG
        AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(daughterCandidate->TrueObject);
        if (trueDaughter && trueDaughter->PDG == 2212){
          nReconProtons++;
        }
        if(trueDaughter && trueDaughter->PDG == 211){
          nRecoPiPlus++;
        }
        if(trueDaughter && trueDaughter->PDG == -211){
          nRecoPiMinus++;
        }

        // Count by chi2/ndf for proton hypothesis
        std::pair<double, int> chi2Proton = pdAnaUtils::Chi2PID(*daughterCandidate, 2212); // Proton PDG
        double chi2ndfProton = (chi2Proton.second > 0) ? chi2Proton.first / chi2Proton.second : 9999.0;
        if(chi2ndfProton < 60.0){
          nReconProtChi2++;
        }
      }
    }
    output.FillVectorVar(k0brothreconprot, (Int_t)nReconProtons);
    output.FillVectorVar(k0brothreconpiplus, (Int_t)nRecoPiPlus);
    output.FillVectorVar(k0brothreconpiminus, (Int_t)nRecoPiMinus);
    output.FillVectorVar(k0brothreconprotchi2, (Int_t)nReconProtChi2);
    if(trueParentCandidate){
      output.FillVectorVar(k0truenbrother, (Int_t)trueParentCandidate->Daughters.size());
    }
  } else {
    // If no parent candidate, fill with default values
    output.FillVectorVar(k0reconbrother, -999);
    output.FillVectorVar(k0brothreconprot, -999);
    output.FillVectorVar(k0brothreconpiplus, -999);
    output.FillVectorVar(k0brothreconpiminus, -999);
    output.FillVectorVar(k0brothreconprotchi2, -999);
    output.FillVectorVar(k0truenbrother, -999);
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtxDaughter1(OutputManager& output, AnaParticlePD* daughterCandidate, AnaVertexPD* vertex){
    // Fill all variables for a single K0 daughter1
  if(daughterCandidate){
    Float_t k0dau1recostartpos_val[3] = {-999.0, -999.0, -999.0};
    k0dau1recostartpos_val[0] = daughterCandidate->PositionStart[0];
    k0dau1recostartpos_val[1] = daughterCandidate->PositionStart[1];
    k0dau1recostartpos_val[2] = daughterCandidate->PositionStart[2];
    output.FillMatrixVarFromArray(k0dau1recostartpos, k0dau1recostartpos_val, 3);

    Float_t k0dau1recoendpos_val[3] = {-999.0, -999.0, -999.0};
    k0dau1recoendpos_val[0] = daughterCandidate->PositionEnd[0];
    k0dau1recoendpos_val[1] = daughterCandidate->PositionEnd[1];
    k0dau1recoendpos_val[2] = daughterCandidate->PositionEnd[2];
    output.FillMatrixVarFromArray(k0dau1recoendpos, k0dau1recoendpos_val, 3);

    Float_t k0dau1recostartdir_val[3] = {-999.0, -999.0, -999.0};
    k0dau1recostartdir_val[0] = daughterCandidate->DirectionStart[0];
    k0dau1recostartdir_val[1] = daughterCandidate->DirectionStart[1];
    k0dau1recostartdir_val[2] = daughterCandidate->DirectionStart[2];
    output.FillMatrixVarFromArray(k0dau1recostartdir, k0dau1recostartdir_val, 3);

    // K0 daughter1 fitted direction from vertex fit
    Float_t k0dau1fitdir_val[3] = {-999.0, -999.0, -999.0};
    if(vertex &&
       vertex->FitParticles.size() >= 1 &&
       vertex->FitParticles[0]){
      k0dau1fitdir_val[0] = vertex->FitParticles[0]->FitDirection[0];
      k0dau1fitdir_val[1] = vertex->FitParticles[0]->FitDirection[1];
      k0dau1fitdir_val[2] = vertex->FitParticles[0]->FitDirection[2];
    }
    output.FillMatrixVarFromArray(k0dau1fitdir, k0dau1fitdir_val, 3);

    Float_t k0dau1recoenddir_val[3] = {-999.0, -999.0, -999.0};
    k0dau1recoenddir_val[0] = daughterCandidate->DirectionEnd[0];
    k0dau1recoenddir_val[1] = daughterCandidate->DirectionEnd[1];
    k0dau1recoenddir_val[2] = daughterCandidate->DirectionEnd[2];
    output.FillMatrixVarFromArray(k0dau1recoenddir, k0dau1recoenddir_val, 3);

    Float_t k0dau1recomom_val = -999.0;
    k0dau1recomom_val = daughterCandidate->Momentum;
    output.FillVectorVar(k0dau1recomom, k0dau1recomom_val);

    Float_t k0dau1recolength_val = -999.0;
    k0dau1recolength_val = daughterCandidate->Length;
    output.FillVectorVar(k0dau1recolength, k0dau1recolength_val);

    Int_t k0dau1recondau_val = -999;
    k0dau1recondau_val = daughterCandidate->Daughters.size();
    output.FillVectorVar(k0dau1recondau, k0dau1recondau_val);

    Int_t k0dau1parid_val = -999;
    k0dau1parid_val = daughterCandidate->ParentID;
    output.FillVectorVar(k0dau1parid, k0dau1parid_val);

    Int_t k0dau1istrack_val = -999;
    k0dau1istrack_val = daughterCandidate->Type;
    output.FillVectorVar(k0dau1istrack, k0dau1istrack_val);

    // Calculate chi2 for pion hypothesis
    std::pair<double, int> chi2ResultPion = pdAnaUtils::Chi2PID(*daughterCandidate, 211);
    Float_t k0dau1chi2pion_val = (Float_t)chi2ResultPion.first;
    Int_t k0dau1nptchi2pion_val = chi2ResultPion.second;
    output.FillVectorVar(k0dau1chi2pion, k0dau1chi2pion_val);
    output.FillVectorVar(k0dau1nptchi2pion, k0dau1nptchi2pion_val);

    // Calculate chi2 for proton hypothesis
    std::pair<double, int> chi2ResultProton = pdAnaUtils::Chi2PID(*daughterCandidate, 2212);
    Float_t k0dau1chi2proton_val = (Float_t)chi2ResultProton.first;
    Int_t k0dau1nptchi2proton_val = chi2ResultProton.second;
    output.FillVectorVar(k0dau1chi2proton, k0dau1chi2proton_val);
    output.FillVectorVar(k0dau1nptchi2proton, k0dau1nptchi2proton_val);

    // Calculate chi2 for kaon hypothesis
    std::pair<double, int> chi2ResultKaon = pdAnaUtils::Chi2PID(*daughterCandidate, 321);
    Float_t k0dau1chi2kaon_val = (Float_t)chi2ResultKaon.first;
    Int_t k0dau1nptchi2kaon_val = chi2ResultKaon.second;
    output.FillVectorVar(k0dau1chi2kaon, k0dau1chi2kaon_val);
    output.FillVectorVar(k0dau1nptchi2kaon, k0dau1nptchi2kaon_val);

    // Calculate average dEdx for different residual range thresholds
    Float_t k0dau1avgdedx5cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 5.0);
    output.FillVectorVar(k0dau1avgdedx5cm, k0dau1avgdedx5cm_val);

    Float_t k0dau1avgdedx10cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 10.0);
    output.FillVectorVar(k0dau1avgdedx10cm, k0dau1avgdedx10cm_val);

    Float_t k0dau1avgdedx15cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 15.0);
    output.FillVectorVar(k0dau1avgdedx15cm, k0dau1avgdedx15cm_val);

    Float_t k0dau1avgdedx25cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 25.0);
    output.FillVectorVar(k0dau1avgdedx25cm, k0dau1avgdedx25cm_val);

    Float_t k0dau1avgdedx50cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 50.0);
    output.FillVectorVar(k0dau1avgdedx50cm, k0dau1avgdedx50cm_val);

    // Number of hits in collection plane
    Int_t k0dau1nhits_val = daughterCandidate->Hits[2].size();
    output.FillVectorVar(k0dau1nhits, k0dau1nhits_val);

    AnaTrueParticlePD* trueDaughterCandidate = static_cast<AnaTrueParticlePD*>(daughterCandidate->TrueObject);
    if(trueDaughterCandidate){
      Float_t k0dau1truestartpos_val[3] = {-999.0, -999.0, -999.0};
      k0dau1truestartpos_val[0] = trueDaughterCandidate->Position[0];
      k0dau1truestartpos_val[1] = trueDaughterCandidate->Position[1];
      k0dau1truestartpos_val[2] = trueDaughterCandidate->Position[2];
      output.FillMatrixVarFromArray(k0dau1truestartpos, k0dau1truestartpos_val, 3);

      Float_t k0dau1trueendpos_val[3] = {-999.0, -999.0, -999.0};
      k0dau1trueendpos_val[0] = trueDaughterCandidate->PositionEnd[0];
      k0dau1trueendpos_val[1] = trueDaughterCandidate->PositionEnd[1];
      k0dau1trueendpos_val[2] = trueDaughterCandidate->PositionEnd[2];
      output.FillMatrixVarFromArray(k0dau1trueendpos, k0dau1trueendpos_val, 3);

      Float_t k0dau1truestartdir_val[3] = {-999.0, -999.0, -999.0};
      k0dau1truestartdir_val[0] = trueDaughterCandidate->Direction[0];
      k0dau1truestartdir_val[1] = trueDaughterCandidate->Direction[1];
      k0dau1truestartdir_val[2] = trueDaughterCandidate->Direction[2];
      output.FillMatrixVarFromArray(k0dau1truestartdir, k0dau1truestartdir_val, 3);


      Float_t k0dau1trueenddir_val[3] = {-999.0, -999.0, -999.0};
      k0dau1trueenddir_val[0] = trueDaughterCandidate->DirectionEnd[0];
      k0dau1trueenddir_val[1] = trueDaughterCandidate->DirectionEnd[1];
      k0dau1trueenddir_val[2] = trueDaughterCandidate->DirectionEnd[2];
      output.FillMatrixVarFromArray(k0dau1trueenddir, k0dau1trueenddir_val, 3);

      Float_t k0dau1truestartmom_val = -999.0;
      k0dau1truestartmom_val = trueDaughterCandidate->Momentum;
      output.FillVectorVar(k0dau1truestartmom, k0dau1truestartmom_val);

      Float_t k0dau1trueendmom_val = -999.0;
      k0dau1trueendmom_val = trueDaughterCandidate->MomentumEnd;
      output.FillVectorVar(k0dau1trueendmom, k0dau1trueendmom_val);

      Float_t k0dau1truelength_val = -999.0;
      k0dau1truelength_val = sqrt(pow(k0dau1trueendpos_val[0] - k0dau1truestartpos_val[0], 2) + pow(k0dau1trueendpos_val[1] - k0dau1truestartpos_val[1], 2) + pow(k0dau1trueendpos_val[2] - k0dau1truestartpos_val[2], 2));
      output.FillVectorVar(k0dau1truelength, k0dau1truelength_val);

      Int_t k0dau1truendau_val = -999;
      k0dau1truendau_val = trueDaughterCandidate->Daughters.size();
      output.FillVectorVar(k0dau1truendau, k0dau1truendau_val);

      Int_t k0dau1trueproc_val = -999;
      k0dau1trueproc_val = static_cast<Int_t>(trueDaughterCandidate->ProcessStart);
      output.FillVectorVar(k0dau1trueproc, k0dau1trueproc_val);

      Int_t k0dau1trueendproc_val = -999;
      k0dau1trueendproc_val = static_cast<Int_t>(trueDaughterCandidate->ProcessEnd);
      output.FillVectorVar(k0dau1trueendproc, k0dau1trueendproc_val);

      Int_t k0dau1truepdg_val = -999;
      k0dau1truepdg_val = trueDaughterCandidate->PDG;
      output.FillVectorVar(k0dau1truepdg, k0dau1truepdg_val);
    } else {
      // Fill with defaults if no true object
      Float_t k0dau1truestartpos_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau1trueendpos_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau1truestartdir_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau1trueenddir_default[3] = {-999.0, -999.0, -999.0};
      output.FillMatrixVarFromArray(k0dau1truestartpos, k0dau1truestartpos_default, 3);
      output.FillMatrixVarFromArray(k0dau1trueendpos, k0dau1trueendpos_default, 3);
      output.FillMatrixVarFromArray(k0dau1truestartdir, k0dau1truestartdir_default, 3);
      output.FillMatrixVarFromArray(k0dau1trueenddir, k0dau1trueenddir_default, 3);
      output.FillVectorVar(k0dau1truestartmom, -999.0);
      output.FillVectorVar(k0dau1trueendmom, -999.0);
      output.FillVectorVar(k0dau1truelength, -999.0);
      output.FillVectorVar(k0dau1truendau, -999);
      output.FillVectorVar(k0dau1trueproc, -999);
      output.FillVectorVar(k0dau1trueendproc, -999);
      output.FillVectorVar(k0dau1truepdg, -999);
    }
  } else {
    // Fill with defaults if no daughter candidate
    Float_t default_3d[3] = {-999.0, -999.0, -999.0};
    output.FillMatrixVarFromArray(k0dau1recostartpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1recoendpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1recostartdir, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1recoenddir, default_3d, 3);
    output.FillVectorVar(k0dau1recomom, -999.0);
    output.FillVectorVar(k0dau1recolength, -999.0);
    output.FillVectorVar(k0dau1recondau, -999);
    output.FillVectorVar(k0dau1parid, -999);
    output.FillVectorVar(k0dau1istrack, -999);
    output.FillVectorVar(k0dau1chi2pion, -999.0);
    output.FillVectorVar(k0dau1nptchi2pion, -999);
    output.FillVectorVar(k0dau1chi2proton, -999.0);
    output.FillVectorVar(k0dau1nptchi2proton, -999);
    output.FillVectorVar(k0dau1chi2kaon, -999.0);
    output.FillVectorVar(k0dau1nptchi2kaon, -999);
    output.FillVectorVar(k0dau1avgdedx5cm, -999.0);
    output.FillVectorVar(k0dau1avgdedx10cm, -999.0);
    output.FillVectorVar(k0dau1avgdedx15cm, -999.0);
    output.FillVectorVar(k0dau1avgdedx25cm, -999.0);
    output.FillVectorVar(k0dau1avgdedx50cm, -999.0);
    output.FillVectorVar(k0dau1nhits, -999);
    output.FillMatrixVarFromArray(k0dau1truestartpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1trueendpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1truestartdir, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau1trueenddir, default_3d, 3);
    output.FillVectorVar(k0dau1truestartmom, -999.0);
    output.FillVectorVar(k0dau1trueendmom, -999.0);
    output.FillVectorVar(k0dau1truelength, -999.0);
    output.FillVectorVar(k0dau1truendau, -999);
    output.FillVectorVar(k0dau1trueproc, -999);
    output.FillVectorVar(k0dau1trueendproc, -999);
    output.FillVectorVar(k0dau1truepdg, -999);
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0vtxDaughter2(OutputManager& output, AnaParticlePD* daughterCandidate, AnaVertexPD* vertex){
    // Fill all variables for a single K0 daughter2
  if(daughterCandidate){
    Float_t k0dau2recostartpos_val[3] = {-999.0, -999.0, -999.0};
    k0dau2recostartpos_val[0] = daughterCandidate->PositionStart[0];
    k0dau2recostartpos_val[1] = daughterCandidate->PositionStart[1];
    k0dau2recostartpos_val[2] = daughterCandidate->PositionStart[2];
    output.FillMatrixVarFromArray(k0dau2recostartpos, k0dau2recostartpos_val, 3);

    Float_t k0dau2recoendpos_val[3] = {-999.0, -999.0, -999.0};
    k0dau2recoendpos_val[0] = daughterCandidate->PositionEnd[0];
    k0dau2recoendpos_val[1] = daughterCandidate->PositionEnd[1];
    k0dau2recoendpos_val[2] = daughterCandidate->PositionEnd[2];
    output.FillMatrixVarFromArray(k0dau2recoendpos, k0dau2recoendpos_val, 3);

    Float_t k0dau2recostartdir_val[3] = {-999.0, -999.0, -999.0};
    k0dau2recostartdir_val[0] = daughterCandidate->DirectionStart[0];
    k0dau2recostartdir_val[1] = daughterCandidate->DirectionStart[1];
    k0dau2recostartdir_val[2] = daughterCandidate->DirectionStart[2];
    output.FillMatrixVarFromArray(k0dau2recostartdir, k0dau2recostartdir_val, 3);

    // K0 daughter2 fitted direction from vertex fit
    Float_t k0dau2fitdir_val[3] = {-999.0, -999.0, -999.0};
    if(vertex &&
       vertex->FitParticles.size() >= 2 &&
       vertex->FitParticles[1]){
      k0dau2fitdir_val[0] = vertex->FitParticles[1]->FitDirection[0];
      k0dau2fitdir_val[1] = vertex->FitParticles[1]->FitDirection[1];
      k0dau2fitdir_val[2] = vertex->FitParticles[1]->FitDirection[2];
    }
    output.FillMatrixVarFromArray(k0dau2fitdir, k0dau2fitdir_val, 3);

    Float_t k0dau2recoenddir_val[3] = {-999.0, -999.0, -999.0};
    k0dau2recoenddir_val[0] = daughterCandidate->DirectionEnd[0];
    k0dau2recoenddir_val[1] = daughterCandidate->DirectionEnd[1];
    k0dau2recoenddir_val[2] = daughterCandidate->DirectionEnd[2];
    output.FillMatrixVarFromArray(k0dau2recoenddir, k0dau2recoenddir_val, 3);

    Float_t k0dau2recomom_val = -999.0;
    k0dau2recomom_val = daughterCandidate->Momentum;
    output.FillVectorVar(k0dau2recomom, k0dau2recomom_val);

    Float_t k0dau2recolength_val = -999.0;
    k0dau2recolength_val = daughterCandidate->Length;
    output.FillVectorVar(k0dau2recolength, k0dau2recolength_val);

    Int_t k0dau2recondau_val = -999;
    k0dau2recondau_val = daughterCandidate->Daughters.size();
    output.FillVectorVar(k0dau2recondau, k0dau2recondau_val);

    Int_t k0dau2parid_val = -999;
    k0dau2parid_val = daughterCandidate->ParentID;
    output.FillVectorVar(k0dau2parid, k0dau2parid_val);

    Int_t k0dau2istrack_val = -999;
    k0dau2istrack_val = daughterCandidate->Type;
    output.FillVectorVar(k0dau2istrack, k0dau2istrack_val);

    // Calculate chi2 for pion hypothesis
    std::pair<double, int> chi2ResultPion = pdAnaUtils::Chi2PID(*daughterCandidate, 211);
    Float_t k0dau2chi2pion_val = (Float_t)chi2ResultPion.first;
    Int_t k0dau2nptchi2pion_val = chi2ResultPion.second;
    output.FillVectorVar(k0dau2chi2pion, k0dau2chi2pion_val);
    output.FillVectorVar(k0dau2nptchi2pion, k0dau2nptchi2pion_val);

    // Calculate chi2 for proton hypothesis
    std::pair<double, int> chi2ResultProton = pdAnaUtils::Chi2PID(*daughterCandidate, 2212);
    Float_t k0dau2chi2proton_val = (Float_t)chi2ResultProton.first;
    Int_t k0dau2nptchi2proton_val = chi2ResultProton.second;
    output.FillVectorVar(k0dau2chi2proton, k0dau2chi2proton_val);
    output.FillVectorVar(k0dau2nptchi2proton, k0dau2nptchi2proton_val);

    // Calculate chi2 for kaon hypothesis
    std::pair<double, int> chi2ResultKaon = pdAnaUtils::Chi2PID(*daughterCandidate, 321);
    Float_t k0dau2chi2kaon_val = (Float_t)chi2ResultKaon.first;
    Int_t k0dau2nptchi2kaon_val = chi2ResultKaon.second;
    output.FillVectorVar(k0dau2chi2kaon, k0dau2chi2kaon_val);
    output.FillVectorVar(k0dau2nptchi2kaon, k0dau2nptchi2kaon_val);

    // Calculate average dEdx for different residual range thresholds
    Float_t k0dau2avgdedx5cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 5.0);
    output.FillVectorVar(k0dau2avgdedx5cm, k0dau2avgdedx5cm_val);

    Float_t k0dau2avgdedx10cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 10.0);
    output.FillVectorVar(k0dau2avgdedx10cm, k0dau2avgdedx10cm_val);

    Float_t k0dau2avgdedx15cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 15.0);
    output.FillVectorVar(k0dau2avgdedx15cm, k0dau2avgdedx15cm_val);

    Float_t k0dau2avgdedx25cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 25.0);
    output.FillVectorVar(k0dau2avgdedx25cm, k0dau2avgdedx25cm_val);

    Float_t k0dau2avgdedx50cm_val = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterCandidate, 50.0);
    output.FillVectorVar(k0dau2avgdedx50cm, k0dau2avgdedx50cm_val);

    // Number of hits in collection plane
    Int_t k0dau2nhits_val = daughterCandidate->Hits[2].size();
    output.FillVectorVar(k0dau2nhits, k0dau2nhits_val);

    AnaTrueParticlePD* trueDaughterCandidate = static_cast<AnaTrueParticlePD*>(daughterCandidate->TrueObject);
    if(trueDaughterCandidate){
      Float_t k0dau2truestartpos_val[3] = {-999.0, -999.0, -999.0};
      k0dau2truestartpos_val[0] = trueDaughterCandidate->Position[0];
      k0dau2truestartpos_val[1] = trueDaughterCandidate->Position[1];
      k0dau2truestartpos_val[2] = trueDaughterCandidate->Position[2];
      output.FillMatrixVarFromArray(k0dau2truestartpos, k0dau2truestartpos_val, 3);

      Float_t k0dau2trueendpos_val[3] = {-999.0, -999.0, -999.0};
      k0dau2trueendpos_val[0] = trueDaughterCandidate->PositionEnd[0];
      k0dau2trueendpos_val[1] = trueDaughterCandidate->PositionEnd[1];
      k0dau2trueendpos_val[2] = trueDaughterCandidate->PositionEnd[2];
      output.FillMatrixVarFromArray(k0dau2trueendpos, k0dau2trueendpos_val, 3);

      Float_t k0dau2truestartdir_val[3] = {-999.0, -999.0, -999.0};
      k0dau2truestartdir_val[0] = trueDaughterCandidate->Direction[0];
      k0dau2truestartdir_val[1] = trueDaughterCandidate->Direction[1];
      k0dau2truestartdir_val[2] = trueDaughterCandidate->Direction[2];
      output.FillMatrixVarFromArray(k0dau2truestartdir, k0dau2truestartdir_val, 3);

      Float_t k0dau2trueenddir_val[3] = {-999.0, -999.0, -999.0};
      k0dau2trueenddir_val[0] = trueDaughterCandidate->DirectionEnd[0];
      k0dau2trueenddir_val[1] = trueDaughterCandidate->DirectionEnd[1];
      k0dau2trueenddir_val[2] = trueDaughterCandidate->DirectionEnd[2];
      output.FillMatrixVarFromArray(k0dau2trueenddir, k0dau2trueenddir_val, 3);

      Float_t k0dau2truestartmom_val = -999.0;
      k0dau2truestartmom_val = trueDaughterCandidate->Momentum;
      output.FillVectorVar(k0dau2truestartmom, k0dau2truestartmom_val);

      Float_t k0dau2trueendmom_val = -999.0;
      k0dau2trueendmom_val = trueDaughterCandidate->MomentumEnd;
      output.FillVectorVar(k0dau2trueendmom, k0dau2trueendmom_val);

      Float_t k0dau2truelength_val = -999.0;
      k0dau2truelength_val = sqrt(pow(k0dau2trueendpos_val[0] - k0dau2truestartpos_val[0], 2) + pow(k0dau2trueendpos_val[1] - k0dau2truestartpos_val[1], 2) + pow(k0dau2trueendpos_val[2] - k0dau2truestartpos_val[2], 2));
      output.FillVectorVar(k0dau2truelength, k0dau2truelength_val);

      Int_t k0dau2truendau_val = -999;
      k0dau2truendau_val = trueDaughterCandidate->Daughters.size();
      output.FillVectorVar(k0dau2truendau, k0dau2truendau_val);

      Int_t k0dau2trueproc_val = -999;
      k0dau2trueproc_val = static_cast<Int_t>(trueDaughterCandidate->ProcessStart);
      output.FillVectorVar(k0dau2trueproc, k0dau2trueproc_val);

      Int_t k0dau2trueendproc_val = -999;
      k0dau2trueendproc_val = static_cast<Int_t>(trueDaughterCandidate->ProcessEnd);
      output.FillVectorVar(k0dau2trueendproc, k0dau2trueendproc_val);

      Int_t k0dau2truepdg_val = -999;
      k0dau2truepdg_val = trueDaughterCandidate->PDG;
      output.FillVectorVar(k0dau2truepdg, k0dau2truepdg_val);
    } else {
      // Fill with defaults if no true object
      Float_t k0dau2truestartpos_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau2trueendpos_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau2truestartdir_default[3] = {-999.0, -999.0, -999.0};
      Float_t k0dau2trueenddir_default[3] = {-999.0, -999.0, -999.0};
      output.FillMatrixVarFromArray(k0dau2truestartpos, k0dau2truestartpos_default, 3);
      output.FillMatrixVarFromArray(k0dau2trueendpos, k0dau2trueendpos_default, 3);
      output.FillMatrixVarFromArray(k0dau2truestartdir, k0dau2truestartdir_default, 3);
      output.FillMatrixVarFromArray(k0dau2trueenddir, k0dau2trueenddir_default, 3);
      output.FillVectorVar(k0dau2truestartmom, -999.0);
      output.FillVectorVar(k0dau2trueendmom, -999.0);
      output.FillVectorVar(k0dau2truelength, -999.0);
      output.FillVectorVar(k0dau2truendau, -999);
      output.FillVectorVar(k0dau2trueproc, -999);
      output.FillVectorVar(k0dau2trueendproc, -999);
      output.FillVectorVar(k0dau2truepdg, -999);
    }
  } else {
    // Fill with defaults if no daughter candidate
    Float_t default_3d[3] = {-999.0, -999.0, -999.0};
    output.FillMatrixVarFromArray(k0dau2recostartpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2recoendpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2recostartdir, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2recoenddir, default_3d, 3);
    output.FillVectorVar(k0dau2recomom, -999.0);
    output.FillVectorVar(k0dau2recolength, -999.0);
    output.FillVectorVar(k0dau2recondau, -999);
    output.FillVectorVar(k0dau2parid, -999);
    output.FillVectorVar(k0dau2istrack, -999);
    output.FillVectorVar(k0dau2chi2pion, -999.0);
    output.FillVectorVar(k0dau2nptchi2pion, -999);
    output.FillVectorVar(k0dau2chi2proton, -999.0);
    output.FillVectorVar(k0dau2nptchi2proton, -999);
    output.FillVectorVar(k0dau2chi2kaon, -999.0);
    output.FillVectorVar(k0dau2nptchi2kaon, -999);
    output.FillVectorVar(k0dau2avgdedx5cm, -999.0);
    output.FillVectorVar(k0dau2avgdedx10cm, -999.0);
    output.FillVectorVar(k0dau2avgdedx15cm, -999.0);
    output.FillVectorVar(k0dau2avgdedx25cm, -999.0);
    output.FillVectorVar(k0dau2avgdedx50cm, -999.0);
    output.FillVectorVar(k0dau2nhits, -999);
    output.FillMatrixVarFromArray(k0dau2truestartpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2trueendpos, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2truestartdir, default_3d, 3);
    output.FillMatrixVarFromArray(k0dau2trueenddir, default_3d, 3);
    output.FillVectorVar(k0dau2truestartmom, -999.0);
    output.FillVectorVar(k0dau2trueendmom, -999.0);
    output.FillVectorVar(k0dau2truelength, -999.0);
    output.FillVectorVar(k0dau2truendau, -999);
    output.FillVectorVar(k0dau2trueproc, -999);
    output.FillVectorVar(k0dau2trueendproc, -999);
    output.FillVectorVar(k0dau2truepdg, -999);
  }
}