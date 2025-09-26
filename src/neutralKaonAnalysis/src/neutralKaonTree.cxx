#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdDataClasses.hxx"
#include "TVector3.h"
#include <cmath>
#include <set>
#include <iostream>

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_Candidates(OutputManager& output, UInt_t nmax){
    //********************************************************************

  // Number of neutral particle candidates
  AddVarI(output, nk0, "number of neutral particle candidates");

  // Unique IDs of neutral particle candidates
  AddVarMaxSizeVI(output, k0id, "unique ID of neutral particle candidates", nk0, nmax);

  // K0 length (distance between parent end position and vertex position)
  AddVarMaxSizeVD(output, k0length, "K0 length - distance between parent end position and vertex position", nk0, nmax);

  // K0 reconstructed end position (vertex position)
  AddVarMaxSizeVD(output, k0recoendpos, "K0 reconstructed end position - vertex position", nk0, nmax);

  // K0 true end position (true start position of one of the vertex particles)
  AddVarMaxSizeVD(output, k0trueendpos, "K0 true end position - true start position of vertex particle", nk0, nmax);

  // K0 daughter true PDG codes (matrix of dimension 2)
  AddVarMI(output, k0dautruepdg, "K0 daughter true PDG codes", nk0, nmax, 2);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_Candidates(OutputManager& output, const std::vector<AnaNeutralParticlePD*>& candidates, const AnaEventB& event){
  //********************************************************************

  // Fill number of candidates
  output.FillVar(nk0, (Int_t)candidates.size());
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_SingleCandidate(OutputManager& output, AnaNeutralParticlePD* candidate){
//********************************************************************

  if (!candidate) return;

  // Fill unique ID of single candidate
  output.FillVectorVar(k0id, candidate->UniqueID);

  // Calculate K0 length (distance between parent end position and vertex position)
  double k0length_val = -999.0;
  if (candidate->Parent && candidate->Vertex) {
    TVector3 parentEnd(candidate->Parent->PositionEnd[0], candidate->Parent->PositionEnd[1], candidate->Parent->PositionEnd[2]);
    TVector3 vertexPos(candidate->Vertex->Position[0], candidate->Vertex->Position[1], candidate->Vertex->Position[2]);
    k0length_val = (parentEnd - vertexPos).Mag();
  }
  output.FillVectorVar(k0length, k0length_val);

  // Fill K0 reconstructed end position (vertex position)
  double k0recoendpos_val = -999.0;
  if (candidate->Vertex) {
    k0recoendpos_val = candidate->Vertex->Position[0]; // X coordinate
  }
  output.FillVectorVar(k0recoendpos, k0recoendpos_val);

  // Fill K0 true end position (true start position of one of the vertex particles)
  double k0trueendpos_val = -999.0;
  if (candidate->Vertex && !candidate->Vertex->Particles.empty()) {
    // Use the true start position of the first particle in the vertex
    AnaParticlePD* firstParticle = candidate->Vertex->Particles[0];
    if (firstParticle && firstParticle->TrueObject) {
      AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(firstParticle->TrueObject);
      if (trueParticle) {
        k0trueendpos_val = trueParticle->Position[0]; // X coordinate (start position)
      }
    }
  }
  output.FillVectorVar(k0trueendpos, k0trueendpos_val);

  // Fill K0 daughter true PDG codes
  Int_t k0dautruepdg_val[2] = {-999, -999};
  if (candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
    AnaParticlePD* particle1 = candidate->Vertex->Particles[0];
    AnaParticlePD* particle2 = candidate->Vertex->Particles[1];

    if (particle1 && particle1->TrueObject) {
      AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
      if (trueParticle1) {
        k0dautruepdg_val[0] = trueParticle1->PDG;
      }
    }

    if (particle2 && particle2->TrueObject) {
      AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);
      if (trueParticle2) {
        k0dautruepdg_val[1] = trueParticle2->PDG;
      }
    }
  }

  // Fill matrix elements individually
  output.FillMatrixVar(k0dautruepdg, k0dautruepdg_val[0], -1, 0);
  output.FillMatrixVar(k0dautruepdg, k0dautruepdg_val[1], -1, 1);

  // output.IncrementCounter(nk0);
}
