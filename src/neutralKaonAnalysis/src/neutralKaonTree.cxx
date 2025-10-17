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

namespace {
  // Helper function to get particle mass from PDG code (in GeV/c^2)
  Float_t GetParticleMass(Int_t pdg) {
    switch(abs(pdg)) {
      case 211:  return 0.13957;   // pi+/-
      case 111:  return 0.134977;  // pi0
      case 2212: return 0.938272;  // proton
      case 2112: return 0.939565;  // neutron
      case 321:  return 0.493677;  // K+/-
      case 310:  return 0.497611;  // K0s
      case 311:  return 0.497611;  // K0
      case 130:  return 0.497611;  // K0L
      case 11:   return 0.000511;  // e+/-
      case 13:   return 0.105658;  // mu+/-
      case 22:   return 0.0;       // gamma
      default:   return -999.0;    // Unknown
    }
  }
}

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
  AddVarMaxSizeVI(output, k0nprotonincreationvtx, "K0 number of parent daughter protons near neutral start", nk0, nmax);
  AddVarMaxSizeVI(output, k0nparticlesincreationvtx, "K0 total number of parent daughters near neutral start", nk0, nmax);
  output.AddMatrixVar(k0creationvtxchi2proton, "k0creationvtxchi2proton", "F", "K0 creation vtx chi2 proton (5 closest)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0creationvtxdistances, "k0creationvtxdistances", "F", "K0 creation vtx distances (5 closest)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0creationvtxtruepdg, "k0creationvtxtruepdg", "I", "K0 creation vtx true PDG (5 closest)", nk0, "nk0", -nmax, 5);
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
  AddVarMaxSize3MF(output, k0partruestartdir, "K0 parent true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecoenddir, "K0 parent reconstructed end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0partrueenddir, "K0 parent true end direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0parrecostartenddir, "K0 parent reconstructed start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruestartenddir, "K0 parent true start-end direction scalar product", nk0, nmax);
  AddVarMaxSizeVF(output, k0parrecolength, "K0 parent reconstructed length", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruelength, "K0 parent true length", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruestartmom, "K0 parent true start momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0partrueendmom, "K0 parent true end momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0partrueendenergy, "K0 parent true end energy", nk0, nmax);
  AddVarMaxSizeVF(output, k0partruestartenergy, "K0 parent true start energy", nk0, nmax);
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
   // Summary counters (backward compatibility)
   AddVarMaxSizeVI(output, k0reconbrother, "Number of K0 brothers", nk0, nmax);
   AddVarMaxSizeVI(output, k0truenbrother, "Number of K0 brothers", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconprot, "Number of K0 brothers that are reconstructed protons", nk0, nmax);
   AddVarMaxSizeVF(output, k0brothrecoprotenergy, "Energy of K0 brother proton with highest momentum", nk0, nmax);
   AddVarMaxSizeVF(output, k0brothrecoprotmom, "Momentum of K0 brother proton with highest momentum", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconpiplus, "Number of K0 brothers that are reconstructed pi+", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconpiminus, "Number of K0 brothers that are reconstructed pi-", nk0, nmax);
   AddVarMaxSizeVI(output, k0brothreconprotchi2, "Number of K0 brothers compatible with proton (chi2/ndf<60)", nk0, nmax);

   // True brothers - ALL from trueParent->Daughters
   AddVarMaxSizeVI(output, k0ntruebroth, "Number of true brothers", nk0, nmax);
   output.AddMatrixVar(k0truebrothpdg, "k0truebrothpdg", "I", "True brothers PDG", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothmomentum, "k0truebrothmomentum", "F", "True brothers momentum", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothenergy, "k0truebrothenergy", "F", "True brothers energy", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothprocessstart, "k0truebrothprocessstart", "I", "True brothers ProcessStart", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothlength, "k0truebrothlength", "F", "True brothers length", nk0, "nk0", -nmax, 100);
   AddVarMaxSizeVF(output, k0truebrothtotalmom, "Total momentum of true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothtotalenergy, "Total energy of true brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0truebrothtruetotaldir, "Resultant true direction of true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothprotonmaxenergy, "Max proton energy in true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothprotonmaxmomentum, "Max proton momentum in true brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0truebrothprotonmaxdir, "Direction of max proton in true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothalign, "Alignment between K+ and K0+proton (all true)", nk0, nmax);

   // True brothers with reco - subset with reconstructed objects
   AddVarMaxSizeVI(output, k0ntruebrothreco, "Number of true brothers with reco", nk0, nmax);
   output.AddMatrixVar(k0truebrothrecotruepdg, "k0truebrothrecotruepdg", "I", "True brothers with reco - true PDG", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecotruemom, "k0truebrothrecotruemom", "F", "True brothers with reco - true momentum", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecotrueenergy, "k0truebrothrecotrueenergy", "F", "True brothers with reco - true energy", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecotrueprocessstart, "k0truebrothrecotrueprocessstart", "I", "True brothers with reco - true ProcessStart", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecopdg, "k0truebrothrecopdg", "I", "True brothers with reco - reco PDG", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecomom, "k0truebrothrecomom", "F", "True brothers with reco - reco momentum", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecoenergy, "k0truebrothrecoenergy", "F", "True brothers with reco - reco energy", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecolength, "k0truebrothrecolength", "F", "True brothers with reco - reco length", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0truebrothrecochi2prot, "k0truebrothrecochi2prot", "F", "True brothers with reco - chi2/ndf proton", nk0, "nk0", -nmax, 100);
   AddVarMaxSizeVF(output, k0truebrothrecotruetotalmom, "Total true momentum of reco true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecotruetotalenergy, "Total true energy of reco true brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0truebrothrecotruetotaldir, "Resultant true direction of reco true brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecoprotonmaxenergy, "Max true proton energy in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecoprotonmaxmomentum, "Max true proton momentum in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecoprotonmaxrecoenergy, "Max true proton reco energy in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecoprotonmaxrecomom, "Max true proton reco momentum in reco brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0truebrothrecoprotonmaxdir, "True direction of max proton in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0truebrothrecoalign, "Alignment between K+ and K0+proton (true parent/K0 + reco proton)", nk0, nmax);

   // Reco brothers - ALL from Parent->Daughters
   AddVarMaxSizeVI(output, k0nrecobroth, "Number of reco brothers", nk0, nmax);
   output.AddMatrixVar(k0recobrothpdg, "k0recobrothpdg", "I", "Reco brothers PDG", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothmom, "k0recobrothmom", "F", "Reco brothers momentum", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothenergy, "k0recobrothenergy", "F", "Reco brothers energy", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothlength, "k0recobrothlength", "F", "Reco brothers length", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothchi2prot, "k0recobrothchi2prot", "F", "Reco brothers chi2/ndf proton", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothtruepdg, "k0recobrothtruepdg", "I", "Reco brothers true PDG", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothtrueenergy, "k0recobrothtrueenergy", "F", "Reco brothers true energy", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothtruemom, "k0recobrothtruemom", "F", "Reco brothers true momentum", nk0, "nk0", -nmax, 100);
   output.AddMatrixVar(k0recobrothtrueprocessstart, "k0recobrothtrueprocessstart", "I", "Reco brothers true ProcessStart", nk0, "nk0", -nmax, 100);
   AddVarMaxSize3MF(output, k0recobrothtruetotaldir, "Resultant true direction of reco brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0recobrothrecototaldir, "Resultant reco direction of reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothprotonmaxenergy, "Max reco proton energy in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothprotonmaxmomentum, "Max reco proton momentum in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothprotonmaxtrueenergy, "Max reco proton true energy in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothprotonmaxtruemom, "Max reco proton true momentum in reco brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0recobrothprotonmaxdir, "Reco direction of max proton in reco brothers", nk0, nmax);
   AddVarMaxSize3MF(output, k0recobrothprotonmaxtruedir, "True direction of max proton in reco brothers", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothalign, "Alignment between K+ and K0+proton (all reco)", nk0, nmax);
   AddVarMaxSizeVF(output, k0recobrothtruealign, "Alignment between K+ and K0+proton (true from reco particles)", nk0, nmax);
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
  AddVarMaxSize3MF(output, k0vtxfitpos, "K0 vertex fitted position (geometric/TMinuit/Kalman)", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxfitdir, "K0 vertex fitted direction (geometric/TMinuit/Kalman)", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxisjustavg, "K0 vertex Pandora used simple average (1) or line intersection (0)", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdegbefore, "K0 vertex degeneracy before scoring", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxdegafter, "K0 vertex degeneracy after scoring", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxnrecopart, "K0 vertex number of unique reco particles", nk0, nmax);
  output.AddMatrixVar(k0vtxdegdistances, "k0vtxdegdistances", "F", "K0 vertex degeneracy distances (5 minimum)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0vtxisolationdistances, "k0vtxisolationdistances", "F", "K0 vertex isolation distances Pandora (5 minimum)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0vtxisolationdistancesfit, "k0vtxisolationdistancesfit", "F", "K0 vertex isolation distances fitted (5 minimum)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0vtxisolationstartdistances, "k0vtxisolationstartdistances", "F", "K0 vertex isolation start distances (5 minimum)", nk0, "nk0", -nmax, 5);
  AddVarMaxSizeVI(output, k0vtxisolnproton, "K0 vertex number of isolation protons", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxisolnpion, "K0 vertex number of isolation pions", nk0, nmax);
  output.AddMatrixVar(k0vtxisolisproton, "k0vtxisolisproton", "I", "K0 vertex isolation is proton flags (5 closest)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0vtxisolchi2proton, "k0vtxisolchi2proton", "F", "K0 vertex isolation chi2 proton (5 closest)", nk0, "nk0", -nmax, 5);
  output.AddMatrixVar(k0vtxisollength, "k0vtxisollength", "F", "K0 vertex isolation lengths (5 closest)", nk0, "nk0", -nmax, 5);
  AddVarMaxSizeVI(output, k0vtxisolislongest, "K0 vertex any isolation particle longer than both vtx particles", nk0, nmax);
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
      neutralKaonTree::FillNeutralKaonVariables_K0Par(output, candidate, event, beam);
      neutralKaonTree::FillNeutralKaonVariables_K0Brother(output, candidate, candidate->Parent, event);
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

    // K0 number of protons in creation vertex
    output.FillVectorVar(k0nprotonincreationvtx, candidate->NProtonInCreationVtx);

    // K0 total number of particles in creation vertex
    output.FillVectorVar(k0nparticlesincreationvtx, candidate->NParticlesInCreationVtx);

    // Fill creation vertex chi2/ndf under proton hypothesis (up to 5)
    Float_t creationChi2Proton[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)candidate->CreationVtxChi2Proton.size()) {
        creationChi2Proton[i] = candidate->CreationVtxChi2Proton[i];
      } else {
        creationChi2Proton[i] = -999.0;
      }
    }
    output.FillMatrixVarFromArray(k0creationvtxchi2proton, creationChi2Proton, 5);

    // Fill creation vertex distances (up to 5)
    Float_t creationDistances[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)candidate->CreationVtxDistances.size()) {
        creationDistances[i] = candidate->CreationVtxDistances[i];
      } else {
        creationDistances[i] = -999.0;
      }
    }
    output.FillMatrixVarFromArray(k0creationvtxdistances, creationDistances, 5);

    // Fill creation vertex true PDG codes (up to 5)
    Int_t creationTruePDG[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)candidate->CreationVtxTruePDG.size()) {
        creationTruePDG[i] = candidate->CreationVtxTruePDG[i];
      } else {
        creationTruePDG[i] = -999;
      }
    }
    output.FillMatrixVarFromArray(k0creationvtxtruepdg, creationTruePDG, 5);

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
            k0trueendmom_val = trueNeutralParticle->MomentumEnd;
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
          k0trueendmom_val = trueEquivalentNeutralParticle->MomentumEnd;
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
void neutralKaonTree::FillNeutralKaonVariables_K0Par(OutputManager& output, AnaNeutralParticlePD* neutralCandidate, const AnaEventB& event, AnaBeamB* beam){
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
      Float_t k0partrueendenergy_val = -999.0;
      Float_t k0partruestartenergy_val = -999.0;

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
      k0partrueendmom_val = trueParentCandidate->MomentumEnd;
      k0partruendau_val = trueParentCandidate->Daughters.size();
      k0partrueproc_val = static_cast<Int_t>(trueParentCandidate->ProcessStart);
      k0partrueendproc_val = static_cast<Int_t>(trueParentCandidate->ProcessEnd);
      k0partruepdg_val = trueParentCandidate->PDG;
      k0partruegeneration_val = trueParentCandidate->Generation;

      // Calculate K0 true end energy: E = sqrt(m^2 + p^2)
      if(k0partrueendmom_val > 0.0 && k0partruepdg_val != -999){
        Float_t parentMass = GetParticleMass(k0partruepdg_val);
        if(parentMass > 0.0){
          k0partrueendenergy_val = sqrt(parentMass*parentMass + k0partrueendmom_val*k0partrueendmom_val);
        }
      }
      // Calculate K0 true start energy: E = sqrt(m^2 + p^2)
      if(k0partruestartmom_val > 0.0 && k0partruepdg_val != -999){
        Float_t parentMass = GetParticleMass(k0partruepdg_val);
        if(parentMass > 0.0){
          k0partruestartenergy_val = sqrt(parentMass*parentMass + k0partruestartmom_val*k0partruestartmom_val);
        }
      }
  }

  output.FillMatrixVarFromArray(k0partruestartpos, k0partruestartpos_val, 3);
  output.FillMatrixVarFromArray(k0partruestartdir, k0partruestartdir_val, 3);
  output.FillMatrixVarFromArray(k0partrueendpos, k0partrueendpos_val, 3);
  output.FillMatrixVarFromArray(k0partrueenddir, k0partrueenddir_val, 3);
    output.FillVectorVar(k0partruestartenddir, k0partruestartenddir_val);
  output.FillVectorVar(k0partruelength, k0partruelength_val);
  output.FillVectorVar(k0partruestartmom, k0partruestartmom_val);
  output.FillVectorVar(k0partrueendmom, k0partrueendmom_val);


  output.FillVectorVar(k0partrueendenergy, k0partrueendenergy_val);
  output.FillVectorVar(k0partruestartenergy, k0partruestartenergy_val);

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

    // Fill fitted position
    output.FillMatrixVarFromArray(k0vtxfitpos, vertex->PositionFit, 3);

    // Fill fitted direction
    output.FillMatrixVarFromArray(k0vtxfitdir, vertex->DirectionFit, 3);

    // Fill IsJustAverage flag
    output.FillVectorVar(k0vtxisjustavg, vertex->IsJustAverage);

    // Fill degeneracy variables
    output.FillVectorVar(k0vtxdegbefore, vertex->DegeneracyBeforeScoring);
    output.FillVectorVar(k0vtxdegafter, vertex->DegeneracyAfterScoring);
    output.FillVectorVar(k0vtxnrecopart, vertex->NRecoParticles);

    // Fill degeneracy distances (up to 5 minimum distances)
    Float_t degDistances[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->DegeneracyDistances.size()) {
        degDistances[i] = vertex->DegeneracyDistances[i];
      } else {
        degDistances[i] = -999.0; // Fill with dummy value if fewer than 5 distances
      }
    }
    output.FillMatrixVarFromArray(k0vtxdegdistances, degDistances, 5);

    // Fill isolation distances (up to 5 minimum distances)
    Float_t isoDistances[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationDistances.size()) {
        isoDistances[i] = vertex->IsolationDistances[i];
      } else {
        isoDistances[i] = -999.0; // Fill with dummy value if fewer than 5 distances
      }
    }
    output.FillMatrixVarFromArray(k0vtxisolationdistances, isoDistances, 5);

    // Fill isolation distances from fitted tracks (up to 5 minimum distances)
    Float_t isoDistancesFit[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationDistancesFit.size()) {
        isoDistancesFit[i] = vertex->IsolationDistancesFit[i];
      } else {
        isoDistancesFit[i] = -999.0; // Fill with dummy value if fewer than 5 distances
      }
    }
    output.FillMatrixVarFromArray(k0vtxisolationdistancesfit, isoDistancesFit, 5);

    // Fill isolation start distances (up to 5 minimum distances)
    Float_t isoStartDistances[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationStartDistances.size()) {
        isoStartDistances[i] = vertex->IsolationStartDistances[i];
      } else {
        isoStartDistances[i] = -999.0; // Fill with dummy value if fewer than 5 distances
      }
    }
    output.FillMatrixVarFromArray(k0vtxisolationstartdistances, isoStartDistances, 5);

    // Fill isolation proton count
    output.FillVectorVar(k0vtxisolnproton, vertex->IsolationNProton);

    // Fill isolation pion count
    output.FillVectorVar(k0vtxisolnpion, vertex->IsolationNPion);

    // Fill isolation is-proton flags (up to 5)
    Int_t isoIsProton[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationIsProton.size()) {
        isoIsProton[i] = vertex->IsolationIsProton[i];
      } else {
        isoIsProton[i] = -999; // Fill with dummy value if fewer than 5
      }
    }
    output.FillMatrixVarFromArray(k0vtxisolisproton, isoIsProton, 5);

    // Fill isolation chi2/ndf under proton hypothesis (up to 5)
    Float_t isoChi2Proton[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationChi2Proton.size()) {
        isoChi2Proton[i] = vertex->IsolationChi2Proton[i];
      } else {
        isoChi2Proton[i] = -999.0; // Fill with dummy value if fewer than 5
      }
    }
    output.FillMatrixVarFromArray(k0vtxisolchi2proton, isoChi2Proton, 5);

    // Fill isolation particle lengths (up to 5)
    Float_t isoLength[5];
    for (int i = 0; i < 5; i++) {
      if (i < (int)vertex->IsolationLength.size()) {
        isoLength[i] = vertex->IsolationLength[i];
      } else {
        isoLength[i] = -999.0; // Fill with dummy value if fewer than 5
      }
    }
    output.FillMatrixVarFromArray(k0vtxisollength, isoLength, 5);

    // Fill isolation is-longest flag
    output.FillVectorVar(k0vtxisolislongest, vertex->IsolationIsLongest);
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_K0Brother(OutputManager& output, AnaNeutralParticlePD* neutralCandidate, AnaParticlePD* parentCandidate, const AnaEventB& event){
    // Fill all variables for a single K0 brother
  if(parentCandidate){
    output.FillVectorVar(k0reconbrother, (Int_t)parentCandidate->Daughters.size());
    AnaTrueParticlePD* trueParentCandidate = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);

    // Get minimum length parameter for counting brothers
    double minBrotherLength = ND::params().GetParameterD("neutralKaonAnalysis.MinBrotherLength");

    // ========== LOOP 1: TRUE BROTHERS (ALL from trueParent->Daughters) ==========
    Int_t nTrueBrothers = 0;
    Float_t trueBrothersTotalMom = 0.0;
    Float_t trueBrothersTotalEnergy = 0.0;
    TVector3 trueBrothersTotalDir(0, 0, 0);
    Float_t trueBrothersProtonMaxEnergy = -999.0;
    Float_t trueBrothersProtonMaxMom = -999.0;
    TVector3 trueBrothersProtonMaxDir(-999.0, -999.0, -999.0);

    Int_t trueBrothersPDG[100];
    Float_t trueBrothersMom[100];
    Float_t trueBrothersEnergy[100];
    Int_t trueBrothersProcess[100];
    Float_t trueBrothersLength[100];

    // Initialize arrays
    for(int i = 0; i < 100; i++){
      trueBrothersPDG[i] = -999;
      trueBrothersMom[i] = -999.0;
      trueBrothersEnergy[i] = -999.0;
      trueBrothersProcess[i] = -999;
      trueBrothersLength[i] = -999.0;
    }

    if(trueParentCandidate){
      int nAllTrueBrothers = std::min((int)trueParentCandidate->Daughters.size(), 100);
      output.FillVectorVar(k0truenbrother, (Int_t)trueParentCandidate->Daughters.size());

      for(int i = 0; i < nAllTrueBrothers; i++){
        Int_t daughterID = trueParentCandidate->Daughters[i];

        // Find true daughter in event
        AnaTrueParticlePD* trueBrother = nullptr;
        for(int j = 0; j < event.nTrueParticles; j++){
          if(event.TrueParticles[j] && event.TrueParticles[j]->ID == daughterID){
            trueBrother = static_cast<AnaTrueParticlePD*>(event.TrueParticles[j]);
            break;
          }
        }

        if(trueBrother){
          trueBrothersPDG[nTrueBrothers] = trueBrother->PDG;
          trueBrothersMom[nTrueBrothers] = trueBrother->Momentum;
          trueBrothersProcess[nTrueBrothers] = static_cast<Int_t>(trueBrother->ProcessStart);

          // Calculate true length from positions
          Float_t dx = trueBrother->PositionEnd[0] - trueBrother->Position[0];
          Float_t dy = trueBrother->PositionEnd[1] - trueBrother->Position[1];
          Float_t dz = trueBrother->PositionEnd[2] - trueBrother->Position[2];
          trueBrothersLength[nTrueBrothers] = sqrt(dx*dx + dy*dy + dz*dz);

          // Calculate energy
          Float_t mass = GetParticleMass(trueBrother->PDG);
          if(mass > 0.0 && trueBrother->Momentum > 0.0){
            Float_t energy = sqrt(mass*mass + trueBrother->Momentum*trueBrother->Momentum);
            trueBrothersEnergy[nTrueBrothers] = energy;
            trueBrothersTotalEnergy += energy;

            // Track max proton
            if(trueBrother->PDG == 2212 && energy > trueBrothersProtonMaxEnergy){
              trueBrothersProtonMaxEnergy = energy;
              trueBrothersProtonMaxMom = trueBrother->Momentum;
              trueBrothersProtonMaxDir.SetXYZ(trueBrother->Direction[0],
                                              trueBrother->Direction[1],
                                              trueBrother->Direction[2]);
            }
          }

          // Accumulate momentum and momentum vector
          trueBrothersTotalMom += trueBrother->Momentum;
          TVector3 momVec(trueBrother->Direction[0] * trueBrother->Momentum,
                          trueBrother->Direction[1] * trueBrother->Momentum,
                          trueBrother->Direction[2] * trueBrother->Momentum);
          trueBrothersTotalDir += momVec;

          nTrueBrothers++;
        }
      }
    }

    // Calculate k0truebrothalign (alignment between K+ and K0+proton using all true quantities)
    Float_t k0truebrothalign_val = -999.0;
    if(neutralCandidate && neutralCandidate->TrueObject && parentCandidate &&
       parentCandidate->TrueObject && trueBrothersProtonMaxMom > 0) {

      AnaTrueParticlePD* trueK0 = static_cast<AnaTrueParticlePD*>(neutralCandidate->TrueObject);
      AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);

      // Parent momentum vector (K+ end momentum)
      TVector3 parentMomVec(trueParent->DirectionEnd[0] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[1] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[2] * trueParent->MomentumEnd);

      // K0 momentum vector (true start momentum)
      TVector3 k0MomVec(trueK0->Direction[0] * trueK0->Momentum,
                        trueK0->Direction[1] * trueK0->Momentum,
                        trueK0->Direction[2] * trueK0->Momentum);

      // Max proton momentum vector
      TVector3 protonMomVec = trueBrothersProtonMaxDir * trueBrothersProtonMaxMom;

      // Sum of daughters (K0 + proton)
      TVector3 daughterSum = k0MomVec + protonMomVec;

      // Cache magnitudes to avoid recalculation (OPTIMIZED)
      Float_t parentMag = parentMomVec.Mag();
      Float_t daughterMag = daughterSum.Mag();

      // Calculate normalized dot product (cosine of angle)
      if(parentMag > 0 && daughterMag > 0) {
        k0truebrothalign_val = parentMomVec.Dot(daughterSum) / (parentMag * daughterMag);
      }
    }

    // ========== LOOP 2: TRUE BROTHERS WITH RECO ==========
    Int_t nTrueBrothersReco = 0;
    Float_t trueBrothersRecoTrueTotalMom = 0.0;
    Float_t trueBrothersRecoTrueTotalEnergy = 0.0;
    TVector3 trueBrothersRecoTrueTotalDir(0, 0, 0);
    Float_t trueBrothersRecoProtonMaxEnergy = -999.0;
    Float_t trueBrothersRecoProtonMaxMom = -999.0;
    Float_t trueBrothersRecoProtonMaxRecoEnergy = -999.0;
    Float_t trueBrothersRecoProtonMaxRecoMom = -999.0;
    TVector3 trueBrothersRecoProtonMaxDir(-999.0, -999.0, -999.0);

    Int_t trueBrothersRecoTruePDG[100];
    Float_t trueBrothersRecoTrueMom[100];
    Float_t trueBrothersRecoTrueEnergy[100];
    Int_t trueBrothersRecoTrueProcess[100];
    Int_t trueBrothersRecoPDG[100];
    Float_t trueBrothersRecoMom[100];
    Float_t trueBrothersRecoEnergy[100];
    Float_t trueBrothersRecoLength[100];
    Float_t trueBrothersRecoChi2Prot[100];

    // Initialize arrays
    for(int i = 0; i < 100; i++){
      trueBrothersRecoTruePDG[i] = -999;
      trueBrothersRecoTrueMom[i] = -999.0;
      trueBrothersRecoTrueEnergy[i] = -999.0;
      trueBrothersRecoTrueProcess[i] = -999;
      trueBrothersRecoPDG[i] = -999;
      trueBrothersRecoMom[i] = -999.0;
      trueBrothersRecoEnergy[i] = -999.0;
      trueBrothersRecoLength[i] = -999.0;
      trueBrothersRecoChi2Prot[i] = -999.0;
    }

    if(trueParentCandidate){
      int nAllTrueBrothers = std::min((int)trueParentCandidate->Daughters.size(), 100);

      for(int i = 0; i < nAllTrueBrothers; i++){
        Int_t daughterID = trueParentCandidate->Daughters[i];

        // Find true daughter
        AnaTrueParticlePD* trueBrother = nullptr;
        for(int j = 0; j < event.nTrueParticles; j++){
          if(event.TrueParticles[j] && event.TrueParticles[j]->ID == daughterID){
            trueBrother = static_cast<AnaTrueParticlePD*>(event.TrueParticles[j]);
            break;
          }
        }

        if(trueBrother){
          // Check if this true brother is reconstructed
          AnaParticlePD* recoBrother = nullptr;
          for(int k = 0; k < event.nParticles; k++){
            if(event.Particles[k] && event.Particles[k]->TrueObject){
              if(event.Particles[k]->TrueObject->ID == trueBrother->ID){
                recoBrother = static_cast<AnaParticlePD*>(event.Particles[k]);
                break;
              }
            }
          }

          if(recoBrother && nTrueBrothersReco < 100){
            // Fill true info
            trueBrothersRecoTruePDG[nTrueBrothersReco] = trueBrother->PDG;
            trueBrothersRecoTrueMom[nTrueBrothersReco] = trueBrother->Momentum;
            trueBrothersRecoTrueProcess[nTrueBrothersReco] = static_cast<Int_t>(trueBrother->ProcessStart);

            // Calculate true energy
            Float_t trueMass = GetParticleMass(trueBrother->PDG);
            if(trueMass > 0.0 && trueBrother->Momentum > 0.0){
              Float_t trueEnergy = sqrt(trueMass*trueMass + trueBrother->Momentum*trueBrother->Momentum);
              trueBrothersRecoTrueEnergy[nTrueBrothersReco] = trueEnergy;
              trueBrothersRecoTrueTotalEnergy += trueEnergy;

              // Track max proton
              if(trueBrother->PDG == 2212 && trueEnergy > trueBrothersRecoProtonMaxEnergy){
                trueBrothersRecoProtonMaxEnergy = trueEnergy;
                trueBrothersRecoProtonMaxMom = trueBrother->Momentum;
                trueBrothersRecoProtonMaxDir.SetXYZ(trueBrother->Direction[0],
                                                    trueBrother->Direction[1],
                                                    trueBrother->Direction[2]);

                // Calculate reco energy for this proton
                Float_t recoMass = GetParticleMass(trueBrother->PDG);
                if(recoMass > 0.0 && recoBrother->Momentum > 0.0){
                  trueBrothersRecoProtonMaxRecoEnergy = sqrt(recoMass*recoMass + recoBrother->Momentum*recoBrother->Momentum);
                  trueBrothersRecoProtonMaxRecoMom = recoBrother->Momentum;
                }
              }
            }

            // Accumulate momentum vector
            trueBrothersRecoTrueTotalMom += trueBrother->Momentum;
            TVector3 momVec(trueBrother->Direction[0] * trueBrother->Momentum,
                            trueBrother->Direction[1] * trueBrother->Momentum,
                            trueBrother->Direction[2] * trueBrother->Momentum);
            trueBrothersRecoTrueTotalDir += momVec;

            // Fill reco info
            trueBrothersRecoPDG[nTrueBrothersReco] = recoBrother->isPandora; // Pandora ID
            trueBrothersRecoMom[nTrueBrothersReco] = recoBrother->Momentum;
            trueBrothersRecoLength[nTrueBrothersReco] = recoBrother->Length;

            // Calculate reco energy
            Float_t recoMass = GetParticleMass(trueBrother->PDG);
            if(recoMass > 0.0 && recoBrother->Momentum > 0.0){
              trueBrothersRecoEnergy[nTrueBrothersReco] = sqrt(recoMass*recoMass + recoBrother->Momentum*recoBrother->Momentum);
            }

            // Chi2 for proton hypothesis
            std::pair<double, int> chi2Prot = pdAnaUtils::Chi2PID(*recoBrother, 2212);
            trueBrothersRecoChi2Prot[nTrueBrothersReco] = (chi2Prot.second > 0) ? chi2Prot.first / chi2Prot.second : -999.0;

            nTrueBrothersReco++;
          }
        }
      }
    }

    // Calculate k0truebrothrecoalign (alignment using true parent/K0 and reco proton)
    Float_t k0truebrothrecoalign_val = -999.0;
    if(neutralCandidate && neutralCandidate->TrueObject && parentCandidate &&
       parentCandidate->TrueObject && trueBrothersRecoProtonMaxRecoMom > 0) {

      AnaTrueParticlePD* trueK0 = static_cast<AnaTrueParticlePD*>(neutralCandidate->TrueObject);
      AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);

      // Parent momentum vector (K+ end momentum - true)
      TVector3 parentMomVec(trueParent->DirectionEnd[0] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[1] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[2] * trueParent->MomentumEnd);

      // K0 momentum vector (true start momentum)
      TVector3 k0MomVec(trueK0->Direction[0] * trueK0->Momentum,
                        trueK0->Direction[1] * trueK0->Momentum,
                        trueK0->Direction[2] * trueK0->Momentum);

      // Max proton momentum vector (reco momentum with true direction)
      TVector3 protonMomVec = trueBrothersRecoProtonMaxDir * trueBrothersRecoProtonMaxRecoMom;

      // Sum of daughters (K0 + proton)
      TVector3 daughterSum = k0MomVec + protonMomVec;

      // Cache magnitudes to avoid recalculation (OPTIMIZED)
      Float_t parentMag = parentMomVec.Mag();
      Float_t daughterMag = daughterSum.Mag();

      // Calculate normalized dot product (cosine of angle)
      if(parentMag > 0 && daughterMag > 0) {
        k0truebrothrecoalign_val = parentMomVec.Dot(daughterSum) / (parentMag * daughterMag);
      }
    }

    // ========== LOOP 3: RECO BROTHERS (ALL from Parent->Daughters) ==========
    Int_t nRecoBrothers = 0;
    TVector3 recoBrothersTrueTotalDir(0, 0, 0);
    TVector3 recoBrothersRecoTotalDir(0, 0, 0);
    Float_t recoBrothersProtonMaxEnergy = -999.0;
    Float_t recoBrothersProtonMaxMom = -999.0;
    Float_t recoBrothersProtonMaxTrueEnergy = -999.0;
    Float_t recoBrothersProtonMaxTrueMom = -999.0;
    TVector3 recoBrothersProtonMaxDir(-999.0, -999.0, -999.0);
    TVector3 recoBrothersProtonMaxTrueDir(-999.0, -999.0, -999.0);

    Int_t recoBrothersPDG[100];
    Float_t recoBrothersMom[100];
    Float_t recoBrothersEnergy[100];
    Float_t recoBrothersLength[100];
    Float_t recoBrothersChi2Prot[100];
    Float_t recoBrothersDir[100][3];  // Cache reco directions for performance
    Int_t recoBrothersTruePDG[100];
    Float_t recoBrothersTrueEnergy[100];
    Float_t recoBrothersTrueMom[100];
    Int_t recoBrothersTrueProcess[100];

    // Initialize arrays
    for(int i = 0; i < 100; i++){
      recoBrothersPDG[i] = -999;
      recoBrothersMom[i] = -999.0;
      recoBrothersEnergy[i] = -999.0;
      recoBrothersLength[i] = -999.0;
      recoBrothersChi2Prot[i] = -999.0;
      recoBrothersDir[i][0] = -999.0;
      recoBrothersDir[i][1] = -999.0;
      recoBrothersDir[i][2] = -999.0;
      recoBrothersTruePDG[i] = -999;
      recoBrothersTrueEnergy[i] = -999.0;
      recoBrothersTrueMom[i] = -999.0;
      recoBrothersTrueProcess[i] = -999;
    }

    UInt_t nReconProtons = 0;
    UInt_t nRecoPiPlus = 0;
    UInt_t nRecoPiMinus = 0;
    UInt_t nReconProtChi2 = 0;
    Float_t maxProtonMom = -999.0;
    Float_t maxProtonEnergy = -999.0;

    for(size_t i = 0; i < parentCandidate->Daughters.size() && nRecoBrothers < 100; i++){
      AnaParticlePD* recoBrother = static_cast<AnaParticlePD*>(parentCandidate->Daughters[i]);
      if(recoBrother){
        // Fill reco info
        recoBrothersPDG[nRecoBrothers] = (recoBrother->ReconPDG[0] != -999) ? recoBrother->ReconPDG[0] : -999;
        recoBrothersMom[nRecoBrothers] = recoBrother->RangeMomentum[1];  // Proton CSDA range momentum
        recoBrothersLength[nRecoBrothers] = recoBrother->Length;

        // Cache direction for performance (avoid re-accessing later)
        recoBrothersDir[nRecoBrothers][0] = recoBrother->DirectionStart[0];
        recoBrothersDir[nRecoBrothers][1] = recoBrother->DirectionStart[1];
        recoBrothersDir[nRecoBrothers][2] = recoBrother->DirectionStart[2];

        // Chi2 for proton hypothesis
        std::pair<double, int> chi2Prot = pdAnaUtils::Chi2PID(*recoBrother, 2212);
        recoBrothersChi2Prot[nRecoBrothers] = (chi2Prot.second > 0) ? chi2Prot.first / chi2Prot.second : -999.0;

        // Accumulate reco momentum vector
        TVector3 recoMomVec(recoBrother->DirectionStart[0] * recoBrother->Momentum,
                            recoBrother->DirectionStart[1] * recoBrother->Momentum,
                            recoBrother->DirectionStart[2] * recoBrother->Momentum);
        recoBrothersRecoTotalDir += recoMomVec;

        // Get true info
        AnaTrueParticlePD* trueBrother = static_cast<AnaTrueParticlePD*>(recoBrother->TrueObject);
        if(trueBrother){
          recoBrothersTruePDG[nRecoBrothers] = trueBrother->PDG;
          recoBrothersTrueMom[nRecoBrothers] = trueBrother->Momentum;
          recoBrothersTrueProcess[nRecoBrothers] = static_cast<Int_t>(trueBrother->ProcessStart);

          // Calculate true energy
          Float_t trueMass = GetParticleMass(trueBrother->PDG);
          if(trueMass > 0.0 && trueBrother->Momentum > 0.0){
            Float_t trueEnergy = sqrt(trueMass*trueMass + trueBrother->Momentum*trueBrother->Momentum);
            recoBrothersTrueEnergy[nRecoBrothers] = trueEnergy;
          }

          // Accumulate true momentum vector
          TVector3 trueMomVec(trueBrother->Direction[0] * trueBrother->Momentum,
                              trueBrother->Direction[1] * trueBrother->Momentum,
                              trueBrother->Direction[2] * trueBrother->Momentum);
          recoBrothersTrueTotalDir += trueMomVec;

          // Track max proton (by reco energy)
          Float_t recoMass = GetParticleMass(trueBrother->PDG);
          if(recoMass > 0.0 && recoBrother->Momentum > 0.0){
            Float_t recoEnergy = sqrt(recoMass*recoMass + recoBrother->Momentum*recoBrother->Momentum);
            recoBrothersEnergy[nRecoBrothers] = recoEnergy;

            if(trueBrother->PDG == 2212 && recoEnergy > recoBrothersProtonMaxEnergy){
              recoBrothersProtonMaxEnergy = recoEnergy;
              recoBrothersProtonMaxMom = recoBrother->Momentum;
              recoBrothersProtonMaxTrueEnergy = recoBrothersTrueEnergy[nRecoBrothers];
              recoBrothersProtonMaxTrueMom = trueBrother->Momentum;
              recoBrothersProtonMaxDir.SetXYZ(recoBrother->DirectionStart[0],
                                              recoBrother->DirectionStart[1],
                                              recoBrother->DirectionStart[2]);
              recoBrothersProtonMaxTrueDir.SetXYZ(trueBrother->Direction[0],
                                                  trueBrother->Direction[1],
                                                  trueBrother->Direction[2]);
            }
          }

          // Summary counters (for backward compatibility)
          if(recoBrother->Length >= minBrotherLength){
            if(trueBrother->PDG == 2212){
          nReconProtons++;
              if(trueBrother->Momentum > maxProtonMom){
                maxProtonMom = trueBrother->Momentum;
                Float_t protonMass = 0.938272;
                maxProtonEnergy = sqrt(protonMass*protonMass + maxProtonMom*maxProtonMom);
              }
            }
            if(trueBrother->PDG == 211) nRecoPiPlus++;
            if(trueBrother->PDG == -211) nRecoPiMinus++;
            if(recoBrothersChi2Prot[nRecoBrothers] < 60.0) nReconProtChi2++;
          }
        }

        nRecoBrothers++;
      }
    }

    // Calculate k0recobrothalign (alignment using all reco quantities) - OPTIMIZED
    Float_t k0recobrothalign_val = -999.0;
    if(neutralCandidate && parentCandidate && nRecoBrothers > 0) {
      // Find max energetic proton in SINGLE pass
      Float_t maxProtonRecoEnergy = -999.0;
      Float_t maxProtonRecoMom = -999.0;
      TVector3 maxProtonRecoDir(-999.0, -999.0, -999.0);
      Float_t minChi2 = 999999.0;
      int minChi2Index = -1;

      // Single loop through all brothers (OPTIMIZED: was 2 loops before)
      for(int i = 0; i < nRecoBrothers; i++){
        // Try PDG first
        if(recoBrothersPDG[i] == 2212 && recoBrothersEnergy[i] > maxProtonRecoEnergy){
          maxProtonRecoEnergy = recoBrothersEnergy[i];
          maxProtonRecoMom = recoBrothersMom[i];
          maxProtonRecoDir.SetXYZ(recoBrothersDir[i][0],
                                  recoBrothersDir[i][1],
                                  recoBrothersDir[i][2]);
        }
        // Track best chi2 as fallback (computed in parallel)
        if(recoBrothersChi2Prot[i] > 0 && recoBrothersChi2Prot[i] < minChi2){
          minChi2 = recoBrothersChi2Prot[i];
          minChi2Index = i;
        }
      }

      // Use chi2 fallback if no PDG match found
      if(maxProtonRecoEnergy < 0 && minChi2Index >= 0){
        maxProtonRecoEnergy = recoBrothersEnergy[minChi2Index];
        maxProtonRecoMom = recoBrothersMom[minChi2Index];
        maxProtonRecoDir.SetXYZ(recoBrothersDir[minChi2Index][0],
                                recoBrothersDir[minChi2Index][1],
                                recoBrothersDir[minChi2Index][2]);
      }

      if(maxProtonRecoMom > 0){
        // Parent momentum magnitude (use true if available, otherwise reco)
        Float_t parentMom = parentCandidate->TrueObject ?
                            static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject)->MomentumEnd :
                            parentCandidate->Momentum;

        // K0 momentum magnitude (use true if available, otherwise fallback)
        Float_t k0Mom = neutralCandidate->TrueObject ?
                        static_cast<AnaTrueParticlePD*>(neutralCandidate->TrueObject)->Momentum :
                        -999.0;

        if(parentMom > 0 && k0Mom > 0){
          // Parent momentum vector (reco end direction with true end momentum)
          TVector3 parentMomVec(parentCandidate->DirectionEnd[0] * parentMom,
                                parentCandidate->DirectionEnd[1] * parentMom,
                                parentCandidate->DirectionEnd[2] * parentMom);

          // K0 momentum vector (reco start direction with true start momentum)
          TVector3 k0MomVec(neutralCandidate->DirectionStart[0] * k0Mom,
                            neutralCandidate->DirectionStart[1] * k0Mom,
                            neutralCandidate->DirectionStart[2] * k0Mom);

          // Proton momentum vector
          TVector3 protonMomVec = maxProtonRecoDir * maxProtonRecoMom;

          // Sum of daughters
          TVector3 daughterSum = k0MomVec + protonMomVec;

          // Cache magnitudes to avoid recalculation (OPTIMIZED)
          Float_t parentMag = parentMomVec.Mag();
          Float_t daughterMag = daughterSum.Mag();

          if(parentMag > 0 && daughterMag > 0){
            k0recobrothalign_val = parentMomVec.Dot(daughterSum) / (parentMag * daughterMag);
          }
        }
      }
    }

    // Calculate k0recobrothtruealign (alignment using true info from reco particles)
    Float_t k0recobrothtruealign_val = -999.0;
    if(neutralCandidate && neutralCandidate->TrueObject && parentCandidate &&
       parentCandidate->TrueObject && recoBrothersProtonMaxTrueMom > 0){

      AnaTrueParticlePD* trueK0 = static_cast<AnaTrueParticlePD*>(neutralCandidate->TrueObject);
      AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(parentCandidate->TrueObject);

      // Parent momentum vector (true end)
      TVector3 parentMomVec(trueParent->DirectionEnd[0] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[1] * trueParent->MomentumEnd,
                            trueParent->DirectionEnd[2] * trueParent->MomentumEnd);

      // K0 momentum vector (true start)
      TVector3 k0MomVec(trueK0->Direction[0] * trueK0->Momentum,
                        trueK0->Direction[1] * trueK0->Momentum,
                        trueK0->Direction[2] * trueK0->Momentum);

      // Max proton momentum vector (from reco brothers with true info)
      TVector3 protonMomVec = recoBrothersProtonMaxTrueDir * recoBrothersProtonMaxTrueMom;

      // Sum of daughters
      TVector3 daughterSum = k0MomVec + protonMomVec;

      // Cache magnitudes to avoid recalculation (OPTIMIZED)
      Float_t parentMag = parentMomVec.Mag();
      Float_t daughterMag = daughterSum.Mag();

      if(parentMag > 0 && daughterMag > 0){
        k0recobrothtruealign_val = parentMomVec.Dot(daughterSum) / (parentMag * daughterMag);
      }
    }

    // Fill summary counters (backward compatibility)
    output.FillVectorVar(k0brothreconprot, (Int_t)nReconProtons);
    output.FillVectorVar(k0brothrecoprotenergy, maxProtonEnergy);
    output.FillVectorVar(k0brothrecoprotmom, maxProtonMom);
    output.FillVectorVar(k0brothreconpiplus, (Int_t)nRecoPiPlus);
    output.FillVectorVar(k0brothreconpiminus, (Int_t)nRecoPiMinus);
    output.FillVectorVar(k0brothreconprotchi2, (Int_t)nReconProtChi2);

    // Fill true brothers arrays
    output.FillVectorVar(k0ntruebroth, nTrueBrothers);
    output.FillMatrixVarFromArray(k0truebrothpdg, trueBrothersPDG, 100);
    output.FillMatrixVarFromArray(k0truebrothmomentum, trueBrothersMom, 100);
    output.FillMatrixVarFromArray(k0truebrothenergy, trueBrothersEnergy, 100);
    output.FillMatrixVarFromArray(k0truebrothprocessstart, trueBrothersProcess, 100);
    output.FillMatrixVarFromArray(k0truebrothlength, trueBrothersLength, 100);
    output.FillVectorVar(k0truebrothtotalmom, trueBrothersTotalMom);
    output.FillVectorVar(k0truebrothtotalenergy, trueBrothersTotalEnergy);

    Float_t trueBrothersTotalDirNorm[3] = {-999.0, -999.0, -999.0};
    if(nTrueBrothers > 0 && trueBrothersTotalDir.Mag() > 0){
      trueBrothersTotalDirNorm[0] = trueBrothersTotalDir.X() / trueBrothersTotalDir.Mag();
      trueBrothersTotalDirNorm[1] = trueBrothersTotalDir.Y() / trueBrothersTotalDir.Mag();
      trueBrothersTotalDirNorm[2] = trueBrothersTotalDir.Z() / trueBrothersTotalDir.Mag();
    }
    output.FillMatrixVarFromArray(k0truebrothtruetotaldir, trueBrothersTotalDirNorm, 3);
    output.FillVectorVar(k0truebrothprotonmaxenergy, trueBrothersProtonMaxEnergy);
    output.FillVectorVar(k0truebrothprotonmaxmomentum, trueBrothersProtonMaxMom);

    Float_t trueBrothersProtonMaxDirArray[3] = {static_cast<Float_t>(trueBrothersProtonMaxDir.X()),
                                                 static_cast<Float_t>(trueBrothersProtonMaxDir.Y()),
                                                 static_cast<Float_t>(trueBrothersProtonMaxDir.Z())};
    output.FillMatrixVarFromArray(k0truebrothprotonmaxdir, trueBrothersProtonMaxDirArray, 3);
    output.FillVectorVar(k0truebrothalign, k0truebrothalign_val);

    // Fill true brothers with reco arrays
    output.FillVectorVar(k0ntruebrothreco, nTrueBrothersReco);
    output.FillMatrixVarFromArray(k0truebrothrecotruepdg, trueBrothersRecoTruePDG, 100);
    output.FillMatrixVarFromArray(k0truebrothrecotruemom, trueBrothersRecoTrueMom, 100);
    output.FillMatrixVarFromArray(k0truebrothrecotrueenergy, trueBrothersRecoTrueEnergy, 100);
    output.FillMatrixVarFromArray(k0truebrothrecotrueprocessstart, trueBrothersRecoTrueProcess, 100);
    output.FillMatrixVarFromArray(k0truebrothrecopdg, trueBrothersRecoPDG, 100);
    output.FillMatrixVarFromArray(k0truebrothrecomom, trueBrothersRecoMom, 100);
    output.FillMatrixVarFromArray(k0truebrothrecoenergy, trueBrothersRecoEnergy, 100);
    output.FillMatrixVarFromArray(k0truebrothrecolength, trueBrothersRecoLength, 100);
    output.FillMatrixVarFromArray(k0truebrothrecochi2prot, trueBrothersRecoChi2Prot, 100);
    output.FillVectorVar(k0truebrothrecotruetotalmom, trueBrothersRecoTrueTotalMom);
    output.FillVectorVar(k0truebrothrecotruetotalenergy, trueBrothersRecoTrueTotalEnergy);

    Float_t trueBrothersRecoTrueTotalDirNorm[3] = {-999.0, -999.0, -999.0};
    if(nTrueBrothersReco > 0 && trueBrothersRecoTrueTotalDir.Mag() > 0){
      trueBrothersRecoTrueTotalDirNorm[0] = trueBrothersRecoTrueTotalDir.X() / trueBrothersRecoTrueTotalDir.Mag();
      trueBrothersRecoTrueTotalDirNorm[1] = trueBrothersRecoTrueTotalDir.Y() / trueBrothersRecoTrueTotalDir.Mag();
      trueBrothersRecoTrueTotalDirNorm[2] = trueBrothersRecoTrueTotalDir.Z() / trueBrothersRecoTrueTotalDir.Mag();
    }
    output.FillMatrixVarFromArray(k0truebrothrecotruetotaldir, trueBrothersRecoTrueTotalDirNorm, 3);
    output.FillVectorVar(k0truebrothrecoprotonmaxenergy, trueBrothersRecoProtonMaxEnergy);
    output.FillVectorVar(k0truebrothrecoprotonmaxmomentum, trueBrothersRecoProtonMaxMom);
    output.FillVectorVar(k0truebrothrecoprotonmaxrecoenergy, trueBrothersRecoProtonMaxRecoEnergy);
    output.FillVectorVar(k0truebrothrecoprotonmaxrecomom, trueBrothersRecoProtonMaxRecoMom);

    Float_t trueBrothersRecoProtonMaxDirArray[3] = {static_cast<Float_t>(trueBrothersRecoProtonMaxDir.X()),
                                                     static_cast<Float_t>(trueBrothersRecoProtonMaxDir.Y()),
                                                     static_cast<Float_t>(trueBrothersRecoProtonMaxDir.Z())};
    output.FillMatrixVarFromArray(k0truebrothrecoprotonmaxdir, trueBrothersRecoProtonMaxDirArray, 3);
    output.FillVectorVar(k0truebrothrecoalign, k0truebrothrecoalign_val);

    // Fill reco brothers arrays
    output.FillVectorVar(k0nrecobroth, nRecoBrothers);
    output.FillMatrixVarFromArray(k0recobrothpdg, recoBrothersPDG, 100);
    output.FillMatrixVarFromArray(k0recobrothmom, recoBrothersMom, 100);
    output.FillMatrixVarFromArray(k0recobrothenergy, recoBrothersEnergy, 100);
    output.FillMatrixVarFromArray(k0recobrothlength, recoBrothersLength, 100);
    output.FillMatrixVarFromArray(k0recobrothchi2prot, recoBrothersChi2Prot, 100);
    output.FillMatrixVarFromArray(k0recobrothtruepdg, recoBrothersTruePDG, 100);
    output.FillMatrixVarFromArray(k0recobrothtrueenergy, recoBrothersTrueEnergy, 100);
    output.FillMatrixVarFromArray(k0recobrothtruemom, recoBrothersTrueMom, 100);
    output.FillMatrixVarFromArray(k0recobrothtrueprocessstart, recoBrothersTrueProcess, 100);

    Float_t recoBrothersTrueTotalDirNorm[3] = {-999.0, -999.0, -999.0};
    if(nRecoBrothers > 0 && recoBrothersTrueTotalDir.Mag() > 0){
      recoBrothersTrueTotalDirNorm[0] = recoBrothersTrueTotalDir.X() / recoBrothersTrueTotalDir.Mag();
      recoBrothersTrueTotalDirNorm[1] = recoBrothersTrueTotalDir.Y() / recoBrothersTrueTotalDir.Mag();
      recoBrothersTrueTotalDirNorm[2] = recoBrothersTrueTotalDir.Z() / recoBrothersTrueTotalDir.Mag();
    }
    output.FillMatrixVarFromArray(k0recobrothtruetotaldir, recoBrothersTrueTotalDirNorm, 3);

    Float_t recoBrothersRecoTotalDirNorm[3] = {-999.0, -999.0, -999.0};
    if(nRecoBrothers > 0 && recoBrothersRecoTotalDir.Mag() > 0){
      recoBrothersRecoTotalDirNorm[0] = recoBrothersRecoTotalDir.X() / recoBrothersRecoTotalDir.Mag();
      recoBrothersRecoTotalDirNorm[1] = recoBrothersRecoTotalDir.Y() / recoBrothersRecoTotalDir.Mag();
      recoBrothersRecoTotalDirNorm[2] = recoBrothersRecoTotalDir.Z() / recoBrothersRecoTotalDir.Mag();
    }
    output.FillMatrixVarFromArray(k0recobrothrecototaldir, recoBrothersRecoTotalDirNorm, 3);
    output.FillVectorVar(k0recobrothprotonmaxenergy, recoBrothersProtonMaxEnergy);
    output.FillVectorVar(k0recobrothprotonmaxmomentum, recoBrothersProtonMaxMom);
    output.FillVectorVar(k0recobrothprotonmaxtrueenergy, recoBrothersProtonMaxTrueEnergy);
    output.FillVectorVar(k0recobrothprotonmaxtruemom, recoBrothersProtonMaxTrueMom);

    Float_t recoBrothersProtonMaxDirArray[3] = {static_cast<Float_t>(recoBrothersProtonMaxDir.X()),
                                                 static_cast<Float_t>(recoBrothersProtonMaxDir.Y()),
                                                 static_cast<Float_t>(recoBrothersProtonMaxDir.Z())};
    output.FillMatrixVarFromArray(k0recobrothprotonmaxdir, recoBrothersProtonMaxDirArray, 3);

    Float_t recoBrothersProtonMaxTrueDirArray[3] = {static_cast<Float_t>(recoBrothersProtonMaxTrueDir.X()),
                                                     static_cast<Float_t>(recoBrothersProtonMaxTrueDir.Y()),
                                                     static_cast<Float_t>(recoBrothersProtonMaxTrueDir.Z())};
    output.FillMatrixVarFromArray(k0recobrothprotonmaxtruedir, recoBrothersProtonMaxTrueDirArray, 3);
    output.FillVectorVar(k0recobrothalign, k0recobrothalign_val);
    output.FillVectorVar(k0recobrothtruealign, k0recobrothtruealign_val);

  } else {
    // If no parent candidate, fill with default values
    output.FillVectorVar(k0reconbrother, -999);
    output.FillVectorVar(k0truenbrother, -999);
    output.FillVectorVar(k0brothreconprot, -999);
    output.FillVectorVar(k0brothrecoprotenergy, -999.0);
    output.FillVectorVar(k0brothrecoprotmom, -999.0);
    output.FillVectorVar(k0brothreconpiplus, -999);
    output.FillVectorVar(k0brothreconpiminus, -999);
    output.FillVectorVar(k0brothreconprotchi2, -999);

    // Fill new variables with defaults
    output.FillVectorVar(k0ntruebroth, -999);
    output.FillVectorVar(k0truebrothtotalmom, -999.0);
    output.FillVectorVar(k0truebrothtotalenergy, -999.0);
    output.FillVectorVar(k0truebrothprotonmaxenergy, -999.0);
    output.FillVectorVar(k0truebrothprotonmaxmomentum, -999.0);
    output.FillVectorVar(k0truebrothalign, -999.0);

    output.FillVectorVar(k0ntruebrothreco, -999);
    output.FillVectorVar(k0truebrothrecotruetotalmom, -999.0);
    output.FillVectorVar(k0truebrothrecotruetotalenergy, -999.0);
    output.FillVectorVar(k0truebrothrecoprotonmaxenergy, -999.0);
    output.FillVectorVar(k0truebrothrecoprotonmaxmomentum, -999.0);
    output.FillVectorVar(k0truebrothrecoprotonmaxrecoenergy, -999.0);
    output.FillVectorVar(k0truebrothrecoprotonmaxrecomom, -999.0);
    output.FillVectorVar(k0truebrothrecoalign, -999.0);

    output.FillVectorVar(k0nrecobroth, -999);
    output.FillVectorVar(k0recobrothprotonmaxenergy, -999.0);
    output.FillVectorVar(k0recobrothprotonmaxmomentum, -999.0);
    output.FillVectorVar(k0recobrothprotonmaxtrueenergy, -999.0);
    output.FillVectorVar(k0recobrothprotonmaxtruemom, -999.0);
    output.FillVectorVar(k0recobrothalign, -999.0);
    output.FillVectorVar(k0recobrothtruealign, -999.0);
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