#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdDataClasses.hxx"
#include "TVector3.h"
#include "AnalysisUtils.hxx"
#include "ParticleId.hxx"
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


  // K0 daughter true PDG codes (matrix of dimension 2)
  AddVarMI(output, k0dautruepdg, "K0 daughter true PDG codes", nk0, nmax, 2);

  // K0 variables
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
  AddVarMaxSizeVI(output, k0truenbrothers, "K0 true number of brothers", nk0, nmax);
  AddVarMaxSizeVI(output, k0truebrotherspdg, "K0 true brothers PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0trueproc, "K0 true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0trueendproc, "K0 true end process", nk0, nmax);
  AddVarMaxSizeVF(output, k0truerecodist, "K0 true-reco distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0impactparameter, "K0 impact parameter", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxoriginaldistance, "K0 vertex original distance", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxminimumdistance, "K0 vertex minimum distance between fitted lines", nk0, nmax);
  AddVarMaxSizeVI(output, k0hastrueobject, "K0 has true object", nk0, nmax);
  AddVarMaxSizeVI(output, k0truepdg, "K0 true PDG", nk0, nmax);

  // K0 daughter variables
  AddVarMaxSize3MF(output, k0daurecostartpos, "K0 daughter reconstructed start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0dautruestartpos, "K0 daughter true start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0daurecoendpos, "K0 daughter reconstructed end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0dautrueendpos, "K0 daughter true end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0daurecostartdir, "K0 daughter reconstructed start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0dautruestartdir, "K0 daughter true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0daurecoenddir, "K0 daughter reconstructed end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0dautrueenddir, "K0 daughter true end direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0daurecolength, "K0 daughter reconstructed length", nk0, nmax);
  AddVarMaxSizeVF(output, k0dautruelength, "K0 daughter true length", nk0, nmax);
  AddVarMaxSizeVF(output, k0dautruestartmom, "K0 daughter true start momentum", nk0, nmax);
  AddVarMaxSizeVF(output, k0dautrueendmom, "K0 daughter true end momentum", nk0, nmax);
  AddVarMaxSizeVI(output, k0dautruendau, "K0 daughter true number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0dautrueproc, "K0 daughter true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0dautrueendproc, "K0 daughter true end process", nk0, nmax);

  // K0 parent variables
  AddVarMaxSize3MF(output, k0parrecostartpos, "K0 parent reconstructed start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0partruestartpos, "K0 parent true start position", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecoendpos, "K0 parent reconstructed end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0partrueendpos, "K0 parent true end position", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecostartdir, "K0 parent reconstructed start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0partruestartdir, "K0 parent true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0parrecoenddir, "K0 parent reconstructed end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0partrueenddir, "K0 parent true end direction", nk0, nmax);
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

  // Vertex system variables
  AddVarMaxSizeVF(output, k0vtxrecomom, "Vertex system reconstructed momentum magnitude", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtruemom, "Vertex system true momentum magnitude", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoenergy, "Vertex system reconstructed energy", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueenergy, "Vertex system true energy", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecomass, "Vertex system reconstructed mass", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtruemass, "Vertex system true mass", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxrecodir, "Vertex system reconstructed direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxtruedir, "Vertex system true direction", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoopening, "Vertex system reconstructed opening angle", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueopening, "Vertex system true opening angle", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoinvariantmass, "Vertex system reconstructed invariant mass", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueinvariantmass, "Vertex system true invariant mass", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoctau, "Vertex system reconstructed ctau", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtruectau, "Vertex system true ctau", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxrecoangle, "Vertex system reconstructed angle with beam", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxtrueangle, "Vertex system true angle with beam", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxnpotpar, "Number of potential parent particles for vertex", nk0, nmax);

  // Variables about the vertex pions only
  AddVarMaxSizeVF(output, k0vtxpplength, "Vertex pi+ particle length", nk0, nmax);
  AddVarMaxSizeVF(output, k0vtxpmlength, "Vertex pi- particle length", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpptruepdg, "Vertex pi+ particle true PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpmtruepdg, "Vertex pi- particle true PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpptruendau, "Vertex pi+ true number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpmtruendau, "Vertex pi- true number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpprecondau, "Vertex pi+ reco number of daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpmrecondau, "Vertex pi- reco number of daughters", nk0, nmax);
  AddVarMI(output, k0vtxpptruepdgdau, "Vertex pi+ true daughters PDG", nk0, nmax, 200);
  AddVarMI(output, k0vtxpmtruepdgdau, "Vertex pi- true daughters PDG", nk0, nmax, 200);
  AddVarMI(output, k0vtxpprecopdgdau, "Vertex pi+ reco daughters PDG", nk0, nmax, 200);
  AddVarMI(output, k0vtxpmrecopdgdau, "Vertex pi- reco daughters PDG", nk0, nmax, 200);
  AddVarMaxSize3MF(output, k0vtxpptruestartdir, "Vertex pi+ true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpmtruestartdir, "Vertex pi- true start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpptrueenddir, "Vertex pi+ true end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpmtrueenddir, "Vertex pi- true end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpprecostartdir, "Vertex pi+ reco start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpmrecostartdir, "Vertex pi- reco start direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpprecoenddir, "Vertex pi+ reco end direction", nk0, nmax);
  AddVarMaxSize3MF(output, k0vtxpmrecoenddir, "Vertex pi- reco end direction", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpptrueproc, "Vertex pi+ true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpmtrueproc, "Vertex pi- true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpptrueendproc, "Vertex pi+ true end process", nk0, nmax);
  AddVarMaxSizeVI(output, k0vtxpmtrueendproc, "Vertex pi- true end process", nk0, nmax);
  AddVarMI(output, k0vtxtruepptruedaupdg, "Vertex pi+ true daughters PDG", nk0, nmax, 200);
  AddVarMI(output, k0vtxtruepmtruedaupdg, "Vertex pi- true daughters PDG", nk0, nmax, 200);

  // Brother variables (proton daughters of K0 parent)
  AddVarMaxSizeVI(output, k0brothnproton, "Number of proton brothers of K0", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothtruepdg, "Proton brothers true PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothrecopdg, "Proton brothers reco PDG", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothtrueproc, "Proton brothers true process", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothtrueendproc, "Proton brothers true end process", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothtruendau, "Proton brothers true n daughters", nk0, nmax);
  AddVarMaxSizeVI(output, k0brothrecondau, "Proton brothers reco n daughters", nk0, nmax);
  AddVarMaxSizeVF(output, k0brothlength, "Proton brothers length", nk0, nmax);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_Candidates(OutputManager& output, const std::vector<AnaNeutralParticlePD*>& candidates, const AnaEventB& event){
  //********************************************************************

  // Fill number of candidates
  output.FillVar(nk0, (Int_t)candidates.size());

  // Loop through each candidate to fill daughter variables
  for (size_t i = 0; i < candidates.size(); i++) {
    AnaNeutralParticlePD* candidate = candidates[i];
    if (!candidate || !candidate->Vertex) continue;

    // Fill daughter variables for this candidate
    if (candidate->Vertex->Particles.size() >= 2) {
      // Daughter 1
      AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[0]);
      if (daughter1) {
        // Daughter 1 reconstructed variables
        Float_t dau1recostartpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau1recostartdir[3] = {-999.0, -999.0, -999.0};
        Float_t dau1recoendpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau1recoenddir[3] = {-999.0, -999.0, -999.0};
        Float_t dau1recolength = -999.0;

        dau1recostartpos[0] = daughter1->PositionStart[0];
        dau1recostartpos[1] = daughter1->PositionStart[1];
        dau1recostartpos[2] = daughter1->PositionStart[2];
        dau1recostartdir[0] = daughter1->DirectionStart[0];
        dau1recostartdir[1] = daughter1->DirectionStart[1];
        dau1recostartdir[2] = daughter1->DirectionStart[2];
        dau1recoendpos[0] = daughter1->PositionEnd[0];
        dau1recoendpos[1] = daughter1->PositionEnd[1];
        dau1recoendpos[2] = daughter1->PositionEnd[2];
        dau1recoenddir[0] = daughter1->DirectionEnd[0];
        dau1recoenddir[1] = daughter1->DirectionEnd[1];
        dau1recoenddir[2] = daughter1->DirectionEnd[2];
        dau1recolength = daughter1->Length;

        output.FillMatrixVarFromArray(k0daurecostartpos, dau1recostartpos, 3);
        output.FillMatrixVarFromArray(k0daurecostartdir, dau1recostartdir, 3);
        output.FillMatrixVarFromArray(k0daurecoendpos, dau1recoendpos, 3);
        output.FillMatrixVarFromArray(k0daurecoenddir, dau1recoenddir, 3);
        output.FillVectorVar(k0daurecolength, dau1recolength);

        // Daughter 1 true variables
        Float_t dau1truestartpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau1truestartdir[3] = {-999.0, -999.0, -999.0};
        Float_t dau1trueendpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau1trueenddir[3] = {-999.0, -999.0, -999.0};
        Float_t dau1truelength = -999.0;
        Float_t dau1truestartmom = -999.0;
        Float_t dau1trueendmom = -999.0;
        Int_t dau1truendau = -999;
        Int_t dau1trueproc = -999;
        Int_t dau1trueendproc = -999;

        if (daughter1->TrueObject) {
          AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
          if (trueDaughter1) {
            dau1truestartpos[0] = trueDaughter1->Position[0];
            dau1truestartpos[1] = trueDaughter1->Position[1];
            dau1truestartpos[2] = trueDaughter1->Position[2];
            dau1truestartdir[0] = trueDaughter1->Direction[0];
            dau1truestartdir[1] = trueDaughter1->Direction[1];
            dau1truestartdir[2] = trueDaughter1->Direction[2];
            dau1trueendpos[0] = trueDaughter1->PositionEnd[0];
            dau1trueendpos[1] = trueDaughter1->PositionEnd[1];
            dau1trueendpos[2] = trueDaughter1->PositionEnd[2];
            dau1trueenddir[0] = trueDaughter1->DirectionEnd[0];
            dau1trueenddir[1] = trueDaughter1->DirectionEnd[1];
            dau1trueenddir[2] = trueDaughter1->DirectionEnd[2];
            dau1truelength = sqrt(pow(trueDaughter1->PositionEnd[0] - trueDaughter1->Position[0], 2) +
                                 pow(trueDaughter1->PositionEnd[1] - trueDaughter1->Position[1], 2) +
                                 pow(trueDaughter1->PositionEnd[2] - trueDaughter1->Position[2], 2));
            dau1truestartmom = trueDaughter1->Momentum;
            dau1trueendmom = trueDaughter1->MomentumEnd;
            dau1truendau = trueDaughter1->Daughters.size();
            dau1trueproc = static_cast<Int_t>(trueDaughter1->ProcessStart);
            dau1trueendproc = static_cast<Int_t>(trueDaughter1->ProcessEnd);
          }
        }

        output.FillMatrixVarFromArray(k0dautruestartpos, dau1truestartpos, 3);
        output.FillMatrixVarFromArray(k0dautruestartdir, dau1truestartdir, 3);
        output.FillMatrixVarFromArray(k0dautrueendpos, dau1trueendpos, 3);
        output.FillMatrixVarFromArray(k0dautrueenddir, dau1trueenddir, 3);
        output.FillVectorVar(k0dautruelength, dau1truelength);
        output.FillVectorVar(k0dautruestartmom, dau1truestartmom);
        output.FillVectorVar(k0dautrueendmom, dau1trueendmom);
        output.FillVectorVar(k0dautruendau, dau1truendau);
        output.FillVectorVar(k0dautrueproc, dau1trueproc);
        output.FillVectorVar(k0dautrueendproc, dau1trueendproc);

        // Daughter 1 true PDG
        Int_t dau1truepdg = -999;
        if (daughter1->TrueObject) {
          AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
          if (trueDaughter1) {
            dau1truepdg = trueDaughter1->PDG;
          }
        }
        output.FillMatrixVarFromArray(k0dautruepdg, &dau1truepdg, 1);
      }

      // Daughter 2
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[1]);
      if (daughter2) {
        // Daughter 2 reconstructed variables
        Float_t dau2recostartpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau2recostartdir[3] = {-999.0, -999.0, -999.0};
        Float_t dau2recoendpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau2recoenddir[3] = {-999.0, -999.0, -999.0};
        Float_t dau2recolength = -999.0;

        dau2recostartpos[0] = daughter2->PositionStart[0];
        dau2recostartpos[1] = daughter2->PositionStart[1];
        dau2recostartpos[2] = daughter2->PositionStart[2];
        dau2recostartdir[0] = daughter2->DirectionStart[0];
        dau2recostartdir[1] = daughter2->DirectionStart[1];
        dau2recostartdir[2] = daughter2->DirectionStart[2];
        dau2recoendpos[0] = daughter2->PositionEnd[0];
        dau2recoendpos[1] = daughter2->PositionEnd[1];
        dau2recoendpos[2] = daughter2->PositionEnd[2];
        dau2recoenddir[0] = daughter2->DirectionEnd[0];
        dau2recoenddir[1] = daughter2->DirectionEnd[1];
        dau2recoenddir[2] = daughter2->DirectionEnd[2];
        dau2recolength = daughter2->Length;

        output.FillMatrixVarFromArray(k0daurecostartpos, dau2recostartpos, 3);
        output.FillMatrixVarFromArray(k0daurecostartdir, dau2recostartdir, 3);
        output.FillMatrixVarFromArray(k0daurecoendpos, dau2recoendpos, 3);
        output.FillMatrixVarFromArray(k0daurecoenddir, dau2recoenddir, 3);
        output.FillVectorVar(k0daurecolength, dau2recolength);

        // Daughter 2 true variables
        Float_t dau2truestartpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau2truestartdir[3] = {-999.0, -999.0, -999.0};
        Float_t dau2trueendpos[3] = {-999.0, -999.0, -999.0};
        Float_t dau2trueenddir[3] = {-999.0, -999.0, -999.0};
        Float_t dau2truelength = -999.0;
        Float_t dau2truestartmom = -999.0;
        Float_t dau2trueendmom = -999.0;
        Int_t dau2truendau = -999;
        Int_t dau2trueproc = -999;
        Int_t dau2trueendproc = -999;

        if (daughter2->TrueObject) {
          AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
          if (trueDaughter2) {
            dau2truestartpos[0] = trueDaughter2->Position[0];
            dau2truestartpos[1] = trueDaughter2->Position[1];
            dau2truestartpos[2] = trueDaughter2->Position[2];
            dau2truestartdir[0] = trueDaughter2->Direction[0];
            dau2truestartdir[1] = trueDaughter2->Direction[1];
            dau2truestartdir[2] = trueDaughter2->Direction[2];
            dau2trueendpos[0] = trueDaughter2->PositionEnd[0];
            dau2trueendpos[1] = trueDaughter2->PositionEnd[1];
            dau2trueendpos[2] = trueDaughter2->PositionEnd[2];
            dau2trueenddir[0] = trueDaughter2->DirectionEnd[0];
            dau2trueenddir[1] = trueDaughter2->DirectionEnd[1];
            dau2trueenddir[2] = trueDaughter2->DirectionEnd[2];
            dau2truelength = sqrt(pow(trueDaughter2->PositionEnd[0] - trueDaughter2->Position[0], 2) +
                                 pow(trueDaughter2->PositionEnd[1] - trueDaughter2->Position[1], 2) +
                                 pow(trueDaughter2->PositionEnd[2] - trueDaughter2->Position[2], 2));
            dau2truestartmom = trueDaughter2->Momentum;
            dau2trueendmom = trueDaughter2->MomentumEnd;
            dau2truendau = trueDaughter2->Daughters.size();
            dau2trueproc = static_cast<Int_t>(trueDaughter2->ProcessStart);
            dau2trueendproc = static_cast<Int_t>(trueDaughter2->ProcessEnd);
          }
        }

        output.FillMatrixVarFromArray(k0dautruestartpos, dau2truestartpos, 3);
        output.FillMatrixVarFromArray(k0dautruestartdir, dau2truestartdir, 3);
        output.FillMatrixVarFromArray(k0dautrueendpos, dau2trueendpos, 3);
        output.FillMatrixVarFromArray(k0dautrueenddir, dau2trueenddir, 3);
        output.FillVectorVar(k0dautruelength, dau2truelength);
        output.FillVectorVar(k0dautruestartmom, dau2truestartmom);
        output.FillVectorVar(k0dautrueendmom, dau2trueendmom);
        output.FillVectorVar(k0dautruendau, dau2truendau);
        output.FillVectorVar(k0dautrueproc, dau2trueproc);
        output.FillVectorVar(k0dautrueendproc, dau2trueendproc);

        // Daughter 2 true PDG
        Int_t dau2truepdg = -999;
        if (daughter2->TrueObject) {
          AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
          if (trueDaughter2) {
            dau2truepdg = trueDaughter2->PDG;
          }
        }
        output.FillMatrixVarFromArray(k0dautruepdg, &dau2truepdg, 1);
      }
    }
  }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_SingleCandidate(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event){
//********************************************************************

  if (!candidate) return;

  // Fill unique ID of single candidate
  output.FillVectorVar(k0id, candidate->UniqueID);

  // ========================================
  // 1. FILL K0 VARIABLES (no loop needed)
  // ========================================

  // K0 reconstructed start position (from AnaNeutralParticlePD inherited from AnaParticleB)
  Float_t k0recostartpos_val[3] = {-999.0, -999.0, -999.0};
  k0recostartpos_val[0] = candidate->PositionStart[0];
  k0recostartpos_val[1] = candidate->PositionStart[1];
  k0recostartpos_val[2] = candidate->PositionStart[2];
  output.FillMatrixVarFromArray(k0recostartpos, k0recostartpos_val, 3);

  // K0 reconstructed start direction
  Float_t k0recostartdir_val[3] = {-999.0, -999.0, -999.0};
  k0recostartdir_val[0] = candidate->DirectionStart[0];
  k0recostartdir_val[1] = candidate->DirectionStart[1];
  k0recostartdir_val[2] = candidate->DirectionStart[2];
  output.FillMatrixVarFromArray(k0recostartdir, k0recostartdir_val, 3);

  // K0 reconstructed end position (vertex position)
  Float_t k0recoendpos_val[3] = {-999.0, -999.0, -999.0};
  if (candidate->Vertex) {
    k0recoendpos_val[0] = candidate->PositionEnd[0];
    k0recoendpos_val[1] = candidate->PositionEnd[1];
    k0recoendpos_val[2] = candidate->PositionEnd[2];
  }
  output.FillMatrixVarFromArray(k0recoendpos, k0recoendpos_val, 3);

  // K0 reconstructed end direction
  Float_t k0recoenddir_val[3] = {-999.0, -999.0, -999.0};
  k0recoenddir_val[0] = candidate->DirectionEnd[0];
  k0recoenddir_val[1] = candidate->DirectionEnd[1];
  k0recoenddir_val[2] = candidate->DirectionEnd[2];
  output.FillMatrixVarFromArray(k0recoenddir, k0recoenddir_val, 3);

  // Calculate reconstructed scalar product
  Float_t k0recostartenddir_val = k0recostartdir_val[0]*k0recoenddir_val[0] + k0recostartdir_val[1]*k0recoenddir_val[1] + k0recostartdir_val[2]*k0recoenddir_val[2];

  // K0 reconstructed length
  Float_t k0recolength_val = sqrt(pow(candidate->PositionEnd[0] - candidate->PositionStart[0], 2) +
                                  pow(candidate->PositionEnd[1] - candidate->PositionStart[1], 2) +
                                  pow(candidate->PositionEnd[2] - candidate->PositionStart[2], 2));
  output.FillVectorVar(k0recolength, k0recolength_val);

  // K0 true variables (need to find associated true particle)
  Float_t k0truestartpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0truestartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueendpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0truestartenddir_val = -999.0;
  Float_t k0trueendpos_val1[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueendpos_val2[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueenddir_val1[3] = {-999.0, -999.0, -999.0};
  Float_t k0trueenddir_val2[3] = {-999.0, -999.0, -999.0};
  Float_t k0truelength_val = -999.0;
  Float_t k0truestartmom_val = -999.0;
  Float_t k0trueendmom_val = -999.0;
  Int_t k0truendau_val = -999;
  Int_t k0truenbrothers_val = -999;
  Int_t k0truebrotherspdg_val[10] = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
  Int_t k0trueproc_val = -999;
  Int_t k0trueendproc_val = -999;
  Int_t k0trueendproc_val1 = -999;
  Int_t k0trueendproc_val2 = -999;
  Int_t k0hastrueobject_val = -999;
  Int_t k0truepdg_val = -999;

  if(candidate->TrueObject){
    AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(candidate->TrueObject);
    if(trueParticle){
      // Fill true quantities from the true particle
      k0truestartpos_val[0] = trueParticle->Position[0];
      k0truestartpos_val[1] = trueParticle->Position[1];
      k0truestartpos_val[2] = trueParticle->Position[2];
      k0truestartdir_val[0] = trueParticle->Direction[0];
      k0truestartdir_val[1] = trueParticle->Direction[1];
      k0truestartdir_val[2] = trueParticle->Direction[2];
      double k0truestartdir_norm = sqrt(k0truestartdir_val[0]*k0truestartdir_val[0] + k0truestartdir_val[1]*k0truestartdir_val[1] + k0truestartdir_val[2]*k0truestartdir_val[2]);
      k0truestartdir_val[0] = k0truestartdir_val[0] / k0truestartdir_norm;
      k0truestartdir_val[1] = k0truestartdir_val[1] / k0truestartdir_norm;
      k0truestartdir_val[2] = k0truestartdir_val[2] / k0truestartdir_norm;

      k0truestartmom_val = trueParticle->Momentum;
      k0trueendmom_val = trueParticle->MomentumEnd;
      k0truenbrothers_val = trueParticle->Daughters.size();
      k0trueproc_val = static_cast<Int_t>(trueParticle->ProcessStart);
      k0trueendproc_val = static_cast<Int_t>(trueParticle->ProcessEnd);
      k0truendau_val = trueParticle->Daughters.size();

        // Fill k0truebrotherspdg array with PDGs of the parent's daughters (siblings)
        if (candidate->Parent && candidate->Parent->TrueObject) {
          AnaTrueParticlePD* parentTrue = static_cast<AnaTrueParticlePD*>(candidate->Parent->TrueObject);
          if (parentTrue) {
            int brotherIndex = 0;
            for (int i = 0; i < parentTrue->Daughters.size() && brotherIndex < 10; i++) {
              // Skip the neutral kaon itself (it's one of the parent's daughters)
              if (parentTrue->Daughters[i] == trueParticle->ID) {
                continue;
              }
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t j = 0; j < nTrueParts; j++) {
                if (trueParticles[j] && trueParticles[j]->ID == parentTrue->Daughters[i]) {
                  AnaTrueParticlePD* daughterTrue = static_cast<AnaTrueParticlePD*>(trueParticles[j]);
                  if (daughterTrue) {
                    k0truebrotherspdg_val[brotherIndex] = daughterTrue->PDG;
                    brotherIndex++;
                  }
                  break;
                }
              }
            }
          }
        }

      k0hastrueobject_val = 1;
      k0truepdg_val = trueParticle->PDG;

      // Calculate true length from true particle's start and end positions
      k0truelength_val = sqrt(pow(trueParticle->Position[0] - trueParticle->PositionEnd[0], 2) +
                              pow(trueParticle->Position[1] - trueParticle->PositionEnd[1], 2) +
                              pow(trueParticle->Position[2] - trueParticle->PositionEnd[2], 2));

      // Set K0 end position and direction from true particle
      k0trueendpos_val[0] = trueParticle->PositionEnd[0];
      k0trueendpos_val[1] = trueParticle->PositionEnd[1];
      k0trueendpos_val[2] = trueParticle->PositionEnd[2];
      k0trueenddir_val[0] = trueParticle->PositionEnd[0] - trueParticle->Position[0];
      k0trueenddir_val[1] = trueParticle->PositionEnd[1] - trueParticle->Position[1];
      k0trueenddir_val[2] = trueParticle->PositionEnd[2] - trueParticle->Position[2];
      double k0trueenddir_norm = sqrt(k0trueenddir_val[0]*k0trueenddir_val[0] + k0trueenddir_val[1]*k0trueenddir_val[1] + k0trueenddir_val[2]*k0trueenddir_val[2]);
      k0trueenddir_val[0] = k0trueenddir_val[0] / k0trueenddir_norm;
      k0trueenddir_val[1] = k0trueenddir_val[1] / k0trueenddir_norm;
      k0trueenddir_val[2] = k0trueenddir_val[2] / k0trueenddir_norm;

      // Calculate and store the scalar product
      k0truestartenddir_val = k0truestartdir_val[0]*k0trueenddir_val[0] + k0truestartdir_val[1]*k0trueenddir_val[1] + k0truestartdir_val[2]*k0trueenddir_val[2];
    }
    else {
      // trueParticle is nullptr, set default values
      k0hastrueobject_val = 0;
      k0truepdg_val = -999;
    }
  }
  else{
    // Get true particle information from parent
    if (candidate->Parent && candidate->Parent->TrueObject) {
      AnaTrueParticlePD* parentTrue = static_cast<AnaTrueParticlePD*>(candidate->Parent->TrueObject);
      if (parentTrue) {
        k0truestartpos_val[0] = parentTrue->PositionEnd[0];
        k0truestartpos_val[1] = parentTrue->PositionEnd[1];
        k0truestartpos_val[2] = parentTrue->PositionEnd[2];
        k0truestartdir_val[0] = candidate->PositionEnd[0] - candidate->PositionStart[0];
        k0truestartdir_val[1] = candidate->PositionEnd[1] - candidate->PositionStart[1];
        k0truestartdir_val[2] = candidate->PositionEnd[2] - candidate->PositionStart[2];

        double k0truestartdir_norm = sqrt(k0truestartdir_val[0]*k0truestartdir_val[0] + k0truestartdir_val[1]*k0truestartdir_val[1] + k0truestartdir_val[2]*k0truestartdir_val[2]);
        k0truestartdir_val[0] = k0truestartdir_val[0] / k0truestartdir_norm;
        k0truestartdir_val[1] = k0truestartdir_val[1] / k0truestartdir_norm;
        k0truestartdir_val[2] = k0truestartdir_val[2] / k0truestartdir_norm;

        k0truestartmom_val = -999;
        k0trueendmom_val = -999;
        k0truenbrothers_val = parentTrue->Daughters.size();

          // Fill k0truebrotherspdg array with PDGs of the parent's daughters (siblings)
          int brotherIndex = 0;
          for (int i = 0; i < parentTrue->Daughters.size() && brotherIndex < 10; i++) {
            // Skip the neutral kaon itself (it's one of the parent's daughters)
            // In this case, we need to find the neutral kaon's ID from the candidate
            if (candidate->TrueObject) {
              AnaTrueParticlePD* candidateTrue = static_cast<AnaTrueParticlePD*>(candidate->TrueObject);
              if (candidateTrue && parentTrue->Daughters[i] == candidateTrue->ID) {
                continue;
              }
            }
            // Find the true particle corresponding to this daughter ID
            AnaTrueParticleB** trueParticles = event.TrueParticles;
            Int_t nTrueParts = event.nTrueParticles;
            for (Int_t j = 0; j < nTrueParts; j++) {
              if (trueParticles[j] && trueParticles[j]->ID == parentTrue->Daughters[i]) {
                AnaTrueParticlePD* daughterTrue = static_cast<AnaTrueParticlePD*>(trueParticles[j]);
                if (daughterTrue) {
                  k0truebrotherspdg_val[brotherIndex] = daughterTrue->PDG;
                  brotherIndex++;
                }
                break;
              }
            }
          }

        k0trueproc_val = static_cast<Int_t>(parentTrue->ProcessStart);
      }
    }

    // Get true particle information from vertex particles
    if (candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
      // First particle
      if (candidate->Vertex->Particles[0] && candidate->Vertex->Particles[0]->TrueObject) {
        AnaTrueParticlePD* particle1True = static_cast<AnaTrueParticlePD*>(candidate->Vertex->Particles[0]->TrueObject);
        if (particle1True) {
          k0trueendpos_val1[0] = particle1True->Position[0];
          k0trueendpos_val1[1] = particle1True->Position[1];
          k0trueendpos_val1[2] = particle1True->Position[2];
          k0trueenddir_val1[0] = particle1True->Direction[0];
          k0trueenddir_val1[1] = particle1True->Direction[1];
          k0trueenddir_val1[2] = particle1True->Direction[2];
          k0trueendproc_val1 = static_cast<Int_t>(particle1True->ProcessStart);
        }
      }

      // Second particle
      if (candidate->Vertex->Particles[1] && candidate->Vertex->Particles[1]->TrueObject) {
        AnaTrueParticlePD* particle2True = static_cast<AnaTrueParticlePD*>(candidate->Vertex->Particles[1]->TrueObject);
        if (particle2True) {
          k0trueendpos_val2[0] = particle2True->Position[0];
          k0trueendpos_val2[1] = particle2True->Position[1];
          k0trueendpos_val2[2] = particle2True->Position[2];
          k0trueenddir_val2[0] = particle2True->Direction[0];
          k0trueenddir_val2[1] = particle2True->Direction[1];
          k0trueenddir_val2[2] = particle2True->Direction[2];
          k0trueendproc_val2 = static_cast<Int_t>(particle2True->ProcessStart);
        }
      }
    }

    // Calculate true length from parent and vertex particles when no TrueObject
    if (candidate->Parent && candidate->Parent->TrueObject &&
        candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
      AnaTrueParticlePD* parentTrue = static_cast<AnaTrueParticlePD*>(candidate->Parent->TrueObject);
      AnaTrueParticlePD* particle1True = static_cast<AnaTrueParticlePD*>(candidate->Vertex->Particles[0]->TrueObject);
      AnaTrueParticlePD* particle2True = static_cast<AnaTrueParticlePD*>(candidate->Vertex->Particles[1]->TrueObject);

      if (parentTrue && particle1True && particle2True) {
        // Calculate length from parent start position to average of vertex particle positions
        Float_t avgEndPos[3] = {
          static_cast<Float_t>((particle1True->Position[0] + particle2True->Position[0]) / 2.0),
          static_cast<Float_t>((particle1True->Position[1] + particle2True->Position[1]) / 2.0),
          static_cast<Float_t>((particle1True->Position[2] + particle2True->Position[2]) / 2.0)
        };

        k0truelength_val = sqrt(pow(parentTrue->PositionEnd[0] - avgEndPos[0], 2) +
                                pow(parentTrue->PositionEnd[1] - avgEndPos[1], 2) +
                                pow(parentTrue->PositionEnd[2] - avgEndPos[2], 2));
      }
    }

    // Set K0 end position and direction as average of the two particles
    k0trueendpos_val[0] = (k0trueendpos_val1[0] + k0trueendpos_val2[0]) / 2.0;
    k0trueendpos_val[1] = (k0trueendpos_val1[1] + k0trueendpos_val2[1]) / 2.0;
    k0trueendpos_val[2] = (k0trueendpos_val1[2] + k0trueendpos_val2[2]) / 2.0;
    k0trueenddir_val[0] = (k0trueenddir_val1[0] + k0trueenddir_val2[0]);
    k0trueenddir_val[1] = (k0trueenddir_val1[1] + k0trueenddir_val2[1]);
    k0trueenddir_val[2] = (k0trueenddir_val1[2] + k0trueenddir_val2[2]);

    double k0trueenddir_norm = sqrt(k0trueenddir_val[0]*k0trueenddir_val[0] + k0trueenddir_val[1]*k0trueenddir_val[1] + k0trueenddir_val[2]*k0trueenddir_val[2]);
    k0trueenddir_val[0] = k0trueenddir_val[0] / k0trueenddir_norm;
    k0trueenddir_val[1] = k0trueenddir_val[1] / k0trueenddir_norm;
    k0trueenddir_val[2] = k0trueenddir_val[2] / k0trueenddir_norm;

    // Calculate true scalar product for the else block
    k0truestartenddir_val = k0truestartdir_val[0]*k0trueenddir_val[0] + k0truestartdir_val[1]*k0trueenddir_val[1] + k0truestartdir_val[2]*k0trueenddir_val[2];

    // Set K0 number of daughters (typically 2 for K0)
    k0truendau_val = 2; // We have 2 vertex particles

    // Set K0 end process
    if(k0trueendproc_val1 == k0trueendproc_val2){
      k0trueendproc_val = k0trueendproc_val1;
    }
    else{
      k0trueendproc_val = -999;
    }
    k0hastrueobject_val = 0;
    k0truepdg_val = -999;
  }

  output.FillMatrixVarFromArray(k0truestartpos, k0truestartpos_val, 3);
  output.FillMatrixVarFromArray(k0truestartdir, k0truestartdir_val, 3);
  output.FillMatrixVarFromArray(k0trueendpos, k0trueendpos_val, 3);
  output.FillMatrixVarFromArray(k0trueenddir, k0trueenddir_val, 3);
  output.FillVectorVar(k0truestartenddir, k0truestartenddir_val);
  output.FillVectorVar(k0recostartenddir, k0recostartenddir_val);
  output.FillVectorVar(k0truelength, k0truelength_val);
  output.FillVectorVar(k0truestartmom, k0truestartmom_val);
  output.FillVectorVar(k0trueendmom, k0trueendmom_val);
  output.FillVectorVar(k0truendau, k0truendau_val);
  output.FillVectorVar(k0truenbrothers, k0truenbrothers_val);
  // Fill brothers PDG array - fill each element individually
  for (int i = 0; i < 10; i++) {
      output.FillVectorVar(k0truebrotherspdg, k0truebrotherspdg_val[i]);
  }
  output.FillVectorVar(k0trueproc, k0trueproc_val);
  output.FillVectorVar(k0trueendproc, k0trueendproc_val);
  output.FillVectorVar(k0hastrueobject, k0hastrueobject_val);
  output.FillVectorVar(k0truepdg, k0truepdg_val);
  // K0 true-reco distance
  Float_t k0truerecodist_val = sqrt(pow(k0trueendpos_val[0] - k0recoendpos_val[0], 2) +
                                   pow(k0trueendpos_val[1] - k0recoendpos_val[1], 2) +
                                   pow(k0trueendpos_val[2] - k0recoendpos_val[2], 2));
  output.FillVectorVar(k0truerecodist, k0truerecodist_val);

  // K0 impact parameter
  output.FillVectorVar(k0impactparameter, candidate->ImpactParameter);

  // K0 vertex original distance
  // K0 vertex minimum distance between fitted lines
  Float_t k0vtxoriginaldistance_val = -999.0;
  Float_t k0vtxminimumdistance_val = -999.0;
  if (candidate->Vertex) {
    k0vtxoriginaldistance_val = candidate->Vertex->OriginalDistance;
    k0vtxminimumdistance_val = candidate->Vertex->MinimumDistance;
  }
  output.FillVectorVar(k0vtxoriginaldistance, k0vtxoriginaldistance_val);
  output.FillVectorVar(k0vtxminimumdistance, k0vtxminimumdistance_val);

  // ========================================
  // K0 MASS CALCULATION
  // ========================================

  // K0 reconstructed mass calculation
  Float_t k0recomass_val = -999.0;
  Float_t k0truemass_val = -999.0;

  // Helper function to calculate invariant mass from daughter particles
  auto calculateInvariantMass = [](Float_t mom1[3], Float_t mom2[3], Float_t mass1, Float_t mass2) -> Float_t {
    // Calculate energies using E² = p² + m²
    Float_t p1_mag_sq = mom1[0]*mom1[0] + mom1[1]*mom1[1] + mom1[2]*mom1[2];
    Float_t p2_mag_sq = mom2[0]*mom2[0] + mom2[1]*mom2[1] + mom2[2]*mom2[2];

    Float_t E1 = sqrt(p1_mag_sq + mass1*mass1);
    Float_t E2 = sqrt(p2_mag_sq + mass2*mass2);

    // Calculate total momentum vector
    Float_t p_total[3] = {mom1[0] + mom2[0], mom1[1] + mom2[1], mom1[2] + mom2[2]};
    Float_t p_total_mag_sq = p_total[0]*p_total[0] + p_total[1]*p_total[1] + p_total[2]*p_total[2];

    // Invariant mass squared: M² = (E₁ + E₂)² - |p⃗₁ + p⃗₂|²
    Float_t mass_sq = (E1 + E2)*(E1 + E2) - p_total_mag_sq;

    // Return mass (sqrt of mass squared), or -999 if negative
    return (mass_sq >= 0.0) ? sqrt(mass_sq) : -999.0;
  };

  // Calculate reconstructed mass if we have vertex particles
  if (candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[0]);
    AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[1]);

    if (daughter1 && daughter2) {
      // For reconstructed mass, assume daughters are pions (π⁺ and π⁻)
      const Float_t pion_mass = 139.57; // MeV/c²

      // Calculate deposited energy using the existing utility function
      Double_t deposited_energy1 = pdAnaUtils::ComputeDepositedEnergy(daughter1);
      Double_t deposited_energy2 = pdAnaUtils::ComputeDepositedEnergy(daughter2);

      // Calculate kinetic energy using the existing utility function
      Float_t kinetic_energy1 = pdAnaUtils::ComputeKineticEnergy(*daughter1);
      Float_t kinetic_energy2 = pdAnaUtils::ComputeKineticEnergy(*daughter2);

      Float_t mom1_reco[3], mom2_reco[3];
      Float_t mass1_reco = pion_mass, mass2_reco = pion_mass;

      // Method 1: Use kinetic energy if available
      if (kinetic_energy1 > 0 && kinetic_energy2 > 0) {
        // Convert kinetic energy from GeV to MeV and calculate momentum
        Float_t KE1_MeV = kinetic_energy1 * 1000.0; // Convert GeV to MeV
        Float_t KE2_MeV = kinetic_energy2 * 1000.0; // Convert GeV to MeV

        // Calculate momentum magnitude: p = sqrt(E² - m²) where E = KE + m
        Float_t E1 = KE1_MeV + pion_mass;
        Float_t E2 = KE2_MeV + pion_mass;
        Float_t p1_mag = sqrt(E1*E1 - pion_mass*pion_mass);
        Float_t p2_mag = sqrt(E2*E2 - pion_mass*pion_mass);

        // Calculate momentum vectors using direction
        mom1_reco[0] = p1_mag * daughter1->DirectionStart[0];
        mom1_reco[1] = p1_mag * daughter1->DirectionStart[1];
        mom1_reco[2] = p1_mag * daughter1->DirectionStart[2];

        mom2_reco[0] = p2_mag * daughter2->DirectionStart[0];
        mom2_reco[1] = p2_mag * daughter2->DirectionStart[1];
        mom2_reco[2] = p2_mag * daughter2->DirectionStart[2];
      }
      // Method 2: Fallback to deposited energy if kinetic energy not available
      else if (deposited_energy1 > 0 && deposited_energy2 > 0) {
        // Use deposited energy as approximation for kinetic energy
        // This is a rough approximation - deposited energy ≈ kinetic energy for stopping particles
        Float_t KE1_MeV = deposited_energy1;
        Float_t KE2_MeV = deposited_energy2;

        // Calculate momentum magnitude: p = sqrt(E² - m²) where E = KE + m
        Float_t E1 = KE1_MeV + pion_mass;
        Float_t E2 = KE2_MeV + pion_mass;
        Float_t p1_mag = sqrt(E1*E1 - pion_mass*pion_mass);
        Float_t p2_mag = sqrt(E2*E2 - pion_mass*pion_mass);

        // Calculate momentum vectors using direction
        mom1_reco[0] = p1_mag * daughter1->DirectionStart[0];
        mom1_reco[1] = p1_mag * daughter1->DirectionStart[1];
        mom1_reco[2] = p1_mag * daughter1->DirectionStart[2];

        mom2_reco[0] = p2_mag * daughter2->DirectionStart[0];
        mom2_reco[1] = p2_mag * daughter2->DirectionStart[1];
        mom2_reco[2] = p2_mag * daughter2->DirectionStart[2];
      }
      // Method 3: Fallback to range-momentum relation if other methods fail
      else {
        // Use range-momentum relation assuming pions
        Float_t p1_mag = pdAnaUtils::ComputeRangeMomentum(daughter1->Length, 211); // 211 = π⁺, returns GeV/c
        Float_t p2_mag = pdAnaUtils::ComputeRangeMomentum(daughter2->Length, 211); // 211 = π⁺, returns GeV/c

        // Convert from GeV/c to MeV/c
        Float_t p1_mag_MeV = p1_mag * 1000.0;
        Float_t p2_mag_MeV = p2_mag * 1000.0;

        // Calculate momentum vectors using direction
        mom1_reco[0] = p1_mag_MeV * daughter1->DirectionStart[0];
        mom1_reco[1] = p1_mag_MeV * daughter1->DirectionStart[1];
        mom1_reco[2] = p1_mag_MeV * daughter1->DirectionStart[2];

        mom2_reco[0] = p2_mag_MeV * daughter2->DirectionStart[0];
        mom2_reco[1] = p2_mag_MeV * daughter2->DirectionStart[1];
        mom2_reco[2] = p2_mag_MeV * daughter2->DirectionStart[2];
      }

      k0recomass_val = calculateInvariantMass(mom1_reco, mom2_reco, mass1_reco, mass2_reco);
    }
  }

  // Calculate true mass if we have true information
  if (candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[0]);
    AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[1]);


    if (daughter1 && daughter2 && daughter1->TrueObject && daughter2->TrueObject) {
      AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
      AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

      if (trueDaughter1 && trueDaughter2) {
        // Get true masses based on PDG codes (in MeV/c²)
        Float_t mass1_true = 0.0, mass2_true = 0.0;

        // Assign masses based on PDG codes (in MeV/c²)
        if (abs(trueDaughter1->PDG) == 211) mass1_true = 139.57;      // π±
        else if (abs(trueDaughter1->PDG) == 321) mass1_true = 493.68; // K±
        else if (abs(trueDaughter1->PDG) == 2212) mass1_true = 938.27; // p
        else if (abs(trueDaughter1->PDG) == 13) mass1_true = 105.66;   // μ±
        else mass1_true = 139.57; // Default to pion mass

        if (abs(trueDaughter2->PDG) == 211) mass2_true = 139.57;      // π±
        else if (abs(trueDaughter2->PDG) == 321) mass2_true = 493.68; // K±
        else if (abs(trueDaughter2->PDG) == 2212) mass2_true = 938.27; // p
        else if (abs(trueDaughter2->PDG) == 13) mass2_true = 105.66;   // μ±
        else mass2_true = 139.57; // Default to pion mass

        // Convert momentum from GeV/c to MeV/c and create momentum vectors
        Float_t p1_mag_MeV = trueDaughter1->Momentum * 1000.0; // Convert GeV/c to MeV/c
        Float_t p2_mag_MeV = trueDaughter2->Momentum * 1000.0; // Convert GeV/c to MeV/c

        Float_t mom1_true[3] = {p1_mag_MeV * trueDaughter1->Direction[0],
                               p1_mag_MeV * trueDaughter1->Direction[1],
                               p1_mag_MeV * trueDaughter1->Direction[2]};

        Float_t mom2_true[3] = {p2_mag_MeV * trueDaughter2->Direction[0],
                               p2_mag_MeV * trueDaughter2->Direction[1],
                               p2_mag_MeV * trueDaughter2->Direction[2]};

        k0truemass_val = calculateInvariantMass(mom1_true, mom2_true, mass1_true, mass2_true);
      }
    }
  }

  output.FillVectorVar(k0recomass, k0recomass_val);
  output.FillVectorVar(k0truemass, k0truemass_val);

  // ========================================
  // 2. FILL K0 PARENT VARIABLES (no loop needed)
  // ========================================

  // K0 parent reconstructed variables
  Float_t k0parrecostartpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0parrecostartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0parrecoendpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0parrecoenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0parrecolength_val = -999.0;

  if (candidate->Parent) {
    k0parrecostartpos_val[0] = candidate->Parent->PositionStart[0];
    k0parrecostartpos_val[1] = candidate->Parent->PositionStart[1];
    k0parrecostartpos_val[2] = candidate->Parent->PositionStart[2];
    k0parrecostartdir_val[0] = candidate->Parent->DirectionStart[0];
    k0parrecostartdir_val[1] = candidate->Parent->DirectionStart[1];
    k0parrecostartdir_val[2] = candidate->Parent->DirectionStart[2];
    k0parrecoendpos_val[0] = candidate->Parent->PositionEnd[0];
    k0parrecoendpos_val[1] = candidate->Parent->PositionEnd[1];
    k0parrecoendpos_val[2] = candidate->Parent->PositionEnd[2];
    k0parrecoenddir_val[0] = candidate->Parent->DirectionEnd[0];
    k0parrecoenddir_val[1] = candidate->Parent->DirectionEnd[1];
    k0parrecoenddir_val[2] = candidate->Parent->DirectionEnd[2];
    k0parrecolength_val = candidate->Parent->Length;
  }

  output.FillMatrixVarFromArray(k0parrecostartpos, k0parrecostartpos_val, 3);
  output.FillMatrixVarFromArray(k0parrecostartdir, k0parrecostartdir_val, 3);
  output.FillMatrixVarFromArray(k0parrecoendpos, k0parrecoendpos_val, 3);
  output.FillMatrixVarFromArray(k0parrecoenddir, k0parrecoenddir_val, 3);
  output.FillVectorVar(k0parrecolength, k0parrecolength_val);

  // K0 parent true variables
  Float_t k0partruestartpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0partruestartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0partrueendpos_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0partrueenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0partruelength_val = -999.0;
  Float_t k0partruestartmom_val = -999.0;
  Float_t k0partrueendmom_val = -999.0;
  Int_t k0partruendau_val = -999;
  Int_t k0parrecondau_val = -999;
  Int_t k0partruepdg_val = -999;
  Int_t k0parrecopdg_val = -999;
  Int_t k0partruepdgdau_val[50] = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
  Int_t k0parrecopdgdau_val[50] = {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};
  Int_t k0parisbeam_val = 0;
  Int_t k0beaminstpdg_val = -999;
  Int_t k0partrueproc_val = -999;
  Int_t k0partrueendproc_val = -999;

  if (candidate->Parent && candidate->Parent->TrueObject) {
    AnaTrueParticlePD* parentTrueParticle = static_cast<AnaTrueParticlePD*>(candidate->Parent->TrueObject);
    if (parentTrueParticle) {
      k0partruestartpos_val[0] = parentTrueParticle->Position[0];
      k0partruestartpos_val[1] = parentTrueParticle->Position[1];
      k0partruestartpos_val[2] = parentTrueParticle->Position[2];
      k0partruestartdir_val[0] = parentTrueParticle->Direction[0];
      k0partruestartdir_val[1] = parentTrueParticle->Direction[1];
      k0partruestartdir_val[2] = parentTrueParticle->Direction[2];
      k0partrueendpos_val[0] = parentTrueParticle->PositionEnd[0];
      k0partrueendpos_val[1] = parentTrueParticle->PositionEnd[1];
      k0partrueendpos_val[2] = parentTrueParticle->PositionEnd[2];
      k0partrueenddir_val[0] = parentTrueParticle->DirectionEnd[0];
      k0partrueenddir_val[1] = parentTrueParticle->DirectionEnd[1];
      k0partrueenddir_val[2] = parentTrueParticle->DirectionEnd[2];
      k0partruelength_val = sqrt(pow(parentTrueParticle->PositionEnd[0] - parentTrueParticle->Position[0], 2) +
                                pow(parentTrueParticle->PositionEnd[1] - parentTrueParticle->Position[1], 2) +
                                pow(parentTrueParticle->PositionEnd[2] - parentTrueParticle->Position[2], 2));
      k0partruestartmom_val = parentTrueParticle->Momentum;
      k0partrueendmom_val = parentTrueParticle->MomentumEnd;
      k0partruendau_val = parentTrueParticle->Daughters.size();
      k0partrueproc_val = static_cast<Int_t>(parentTrueParticle->ProcessStart);
      k0partrueendproc_val = static_cast<Int_t>(parentTrueParticle->ProcessEnd);
      k0partruepdg_val = parentTrueParticle->PDG;

      // Fill true PDG of parent's daughters
      for (int i = 0; i < parentTrueParticle->Daughters.size() && i < 50; i++) {
        // Find the true particle corresponding to this daughter ID
        AnaTrueParticleB** trueParticles = event.TrueParticles;
        Int_t nTrueParts = event.nTrueParticles;
        for (Int_t j = 0; j < nTrueParts; j++) {
          if (trueParticles[j] && trueParticles[j]->ID == parentTrueParticle->Daughters[i]) {
            AnaTrueParticlePD* daughterTrue = static_cast<AnaTrueParticlePD*>(trueParticles[j]);
            if (daughterTrue) {
              k0partruepdgdau_val[i] = daughterTrue->PDG;
            }
            break;
          }
        }
      }
    }
  }

  // Fill reconstructed parent information
  if (candidate->Parent) {
    k0parrecondau_val = candidate->Parent->Daughters.size();

    // Get reconstructed PDG from ReconPDG (use plane 2 as default, or first non-default value)
    if (candidate->Parent->ReconPDG[2] != -999) {
      k0parrecopdg_val = candidate->Parent->ReconPDG[2];
    } else if (candidate->Parent->ReconPDG[0] != -999) {
      k0parrecopdg_val = candidate->Parent->ReconPDG[0];
    } else if (candidate->Parent->ReconPDG[1] != -999) {
      k0parrecopdg_val = candidate->Parent->ReconPDG[1];
    }

    // Fill reconstructed PDG of parent's daughters
    for (int i = 0; i < candidate->Parent->Daughters.size() && i < 50; i++) {
      AnaParticlePD* daughter = static_cast<AnaParticlePD*>(candidate->Parent->Daughters[i]);
      if (daughter) {
        // Use plane 2 as default, or first non-default value
        if (daughter->ReconPDG[2] != -999) {
          k0parrecopdgdau_val[i] = daughter->ReconPDG[2];
        } else if (daughter->ReconPDG[0] != -999) {
          k0parrecopdgdau_val[i] = daughter->ReconPDG[0];
        } else if (daughter->ReconPDG[1] != -999) {
          k0parrecopdgdau_val[i] = daughter->ReconPDG[1];
        }
      }
    }
  }

  // Check if parent is beam particle and get beam instrumentation PDG
  if (candidate->Parent) {
    // Get the event to access beam information (need const_cast to access Beam)
    const AnaEventPD* eventPD = static_cast<const AnaEventPD*>(&event);
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(eventPD->Beam);
    if (beam && beam->BeamParticle) {
      if (candidate->Parent == beam->BeamParticle) {
        k0parisbeam_val = 1;
      }
      // Get beam instrumentation PDG (first PDG in the PDGs vector if available)
      if (beam->PDGs.size() > 0) {
        k0beaminstpdg_val = beam->PDGs[0];
      }
    }
  }

  output.FillMatrixVarFromArray(k0partruestartpos, k0partruestartpos_val, 3);
  output.FillMatrixVarFromArray(k0partruestartdir, k0partruestartdir_val, 3);
  output.FillMatrixVarFromArray(k0partrueendpos, k0partrueendpos_val, 3);
  output.FillMatrixVarFromArray(k0partrueenddir, k0partrueenddir_val, 3);
  output.FillVectorVar(k0partruelength, k0partruelength_val);
  output.FillVectorVar(k0partruestartmom, k0partruestartmom_val);
  output.FillVectorVar(k0partrueendmom, k0partrueendmom_val);
  output.FillVectorVar(k0partruendau, k0partruendau_val);
  output.FillVectorVar(k0parrecondau, k0parrecondau_val);
  output.FillVectorVar(k0partruepdg, k0partruepdg_val);
  output.FillVectorVar(k0parrecopdg, k0parrecopdg_val);

  // Fill parent daughters PDG arrays
  for (int i = 0; i < 50; i++) {
    output.FillVectorVar(k0partruepdgdau, k0partruepdgdau_val[i]);
  }
  for (int i = 0; i < 50; i++) {
    output.FillVectorVar(k0parrecopdgdau, k0parrecopdgdau_val[i]);
  }

  output.FillVectorVar(k0parisbeam, k0parisbeam_val);
  output.FillVectorVar(k0beaminstpdg, k0beaminstpdg_val);
  output.FillVectorVar(k0partrueproc, k0partrueproc_val);
  output.FillVectorVar(k0partrueendproc, k0partrueendproc_val);

  // ========================================
  // 3. FILL VERTEX SYSTEM VARIABLES
  // ========================================

  // Initialize vertex system variables
  Float_t k0vtxrecomom_val = -999.0;
  Float_t k0vtxtruemom_val = -999.0;
  Float_t k0vtxrecoenergy_val = -999.0;
  Float_t k0vtxtrueenergy_val = -999.0;
  Float_t k0vtxrecomass_val = -999.0;
  Float_t k0vtxtruemass_val = -999.0;
  Float_t k0vtxrecodir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxtruedir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxrecoopening_val = -999.0;
  Float_t k0vtxtrueopening_val = -999.0;
  Float_t k0vtxrecoinvariantmass_val = -999.0;
  Float_t k0vtxtrueinvariantmass_val = -999.0;
  Float_t k0vtxrecoctau_val = -999.0;
  Float_t k0vtxtruectau_val = -999.0;
  Float_t k0vtxrecoangle_val = -999.0;
  Float_t k0vtxtrueangle_val = -999.0;
  Int_t k0vtxnpotpar_val = 0;
  Float_t k0vtxpplength_val = -999.0;
  Float_t k0vtxpmlength_val = -999.0;
  Int_t k0vtxpptruepdg_val = -999;
  Int_t k0vtxpmtruepdg_val = -999;
  Int_t k0vtxpptruendau_val = -999;
  Int_t k0vtxpmtruendau_val = -999;
  Int_t k0vtxpprecondau_val = -999;
  Int_t k0vtxpmrecondau_val = -999;
  Int_t k0vtxpptruepdgdau_val[200];
  Int_t k0vtxpmtruepdgdau_val[200];
  Int_t k0vtxpprecopdgdau_val[200];
  Int_t k0vtxpmrecopdgdau_val[200];

  // Initialize daughter PDG arrays
  for (int i = 0; i < 200; i++) {
    k0vtxpptruepdgdau_val[i] = -999;
    k0vtxpmtruepdgdau_val[i] = -999;
    k0vtxpprecopdgdau_val[i] = -999;
    k0vtxpmrecopdgdau_val[i] = -999;
  }
  Float_t k0vtxpptruestartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpmtruestartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpptrueenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpmtrueenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpprecostartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpmrecostartdir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpprecoenddir_val[3] = {-999.0, -999.0, -999.0};
  Float_t k0vtxpmrecoenddir_val[3] = {-999.0, -999.0, -999.0};
  Int_t k0vtxpptrueproc_val = -999;
  Int_t k0vtxpmtrueproc_val = -999;
  Int_t k0vtxpptrueendproc_val = -999;
  Int_t k0vtxpmtrueendproc_val = -999;
  Int_t k0vtxtruepptruedaupdg_val[200];
  Int_t k0vtxtruepmtruedaupdg_val[200];

  // Brother variables (proton daughters of K0 parent)
  Int_t k0brothnproton_val = 0;
  Int_t k0brothtruepdg_val = -999;
  Int_t k0brothrecopdg_val = -999;
  Int_t k0brothtrueproc_val = -999;
  Int_t k0brothtrueendproc_val = -999;
  Int_t k0brothtruendau_val = -999;
  Int_t k0brothrecondau_val = -999;
  Float_t k0brothlength_val = -999.0;

  // Initialize daughter PDG arrays
  for (int i = 0; i < 200; i++) {
    k0vtxtruepptruedaupdg_val[i] = -999;
    k0vtxtruepmtruedaupdg_val[i] = -999;
  }

  // Calculate vertex system variables if we have vertex particles
  if (candidate->Vertex && candidate->Vertex->Particles.size() >= 2) {
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[0]);
    AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(candidate->Vertex->Particles[1]);

    if (daughter1 && daughter2) {
      // Get the event to access beam and all particles
      const AnaEventPD* eventPD = static_cast<const AnaEventPD*>(&event);

      // ========================================
      // IDENTIFY PI+ AND PI- PARTICLES FROM TRUE INFO
      // ========================================

      // Check true PDG of daughter particles
      Int_t dau1_truepdg = -999;
      Int_t dau2_truepdg = -999;

      if (daughter1->TrueObject) {
        AnaTrueParticlePD* trueDau1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
        if (trueDau1) {
          dau1_truepdg = trueDau1->PDG;
        }
      }

      if (daughter2->TrueObject) {
        AnaTrueParticlePD* trueDau2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
        if (trueDau2) {
          dau2_truepdg = trueDau2->PDG;
        }
      }

      // Identify which particle is pi+ (211) and which is pi- (-211)
      AnaParticlePD* piPlus = nullptr;
      AnaParticlePD* piMinus = nullptr;

      if (dau1_truepdg == 211) {
        piPlus = daughter1;
        k0vtxpptruepdg_val = dau1_truepdg;
        // Get process information for pi+
        if (daughter1->TrueObject) {
          AnaTrueParticlePD* trueDau1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
          if (trueDau1) {
            k0vtxpptrueproc_val = static_cast<Int_t>(trueDau1->ProcessStart);
            k0vtxpptrueendproc_val = static_cast<Int_t>(trueDau1->ProcessEnd);

            // Fill daughter PDGs for pi+
            for (int j = 0; j < trueDau1->Daughters.size() && j < 200; j++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t k = 0; k < nTrueParts; k++) {
                if (trueParticles[k] && trueParticles[k]->ID == trueDau1->Daughters[j]) {
                  AnaTrueParticlePD* trueGrandDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[k]);
                  if (trueGrandDaughter) {
                    k0vtxtruepptruedaupdg_val[j] = trueGrandDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      } else if (dau1_truepdg == -211) {
        piMinus = daughter1;
        k0vtxpmtruepdg_val = dau1_truepdg;
        // Get process information for pi-
        if (daughter1->TrueObject) {
          AnaTrueParticlePD* trueDau1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
          if (trueDau1) {
            k0vtxpmtrueproc_val = static_cast<Int_t>(trueDau1->ProcessStart);
            k0vtxpmtrueendproc_val = static_cast<Int_t>(trueDau1->ProcessEnd);

            // Fill daughter PDGs for pi-
            for (int j = 0; j < trueDau1->Daughters.size() && j < 200; j++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t k = 0; k < nTrueParts; k++) {
                if (trueParticles[k] && trueParticles[k]->ID == trueDau1->Daughters[j]) {
                  AnaTrueParticlePD* trueGrandDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[k]);
                  if (trueGrandDaughter) {
                    k0vtxtruepmtruedaupdg_val[j] = trueGrandDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      if (dau2_truepdg == 211) {
        piPlus = daughter2;
        k0vtxpptruepdg_val = dau2_truepdg;
        // Get process information for pi+
        if (daughter2->TrueObject) {
          AnaTrueParticlePD* trueDau2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
          if (trueDau2) {
            k0vtxpptrueproc_val = static_cast<Int_t>(trueDau2->ProcessStart);
            k0vtxpptrueendproc_val = static_cast<Int_t>(trueDau2->ProcessEnd);

            // Fill daughter PDGs for pi+
            for (int j = 0; j < trueDau2->Daughters.size() && j < 200; j++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t k = 0; k < nTrueParts; k++) {
                if (trueParticles[k] && trueParticles[k]->ID == trueDau2->Daughters[j]) {
                  AnaTrueParticlePD* trueGrandDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[k]);
                  if (trueGrandDaughter) {
                    k0vtxtruepptruedaupdg_val[j] = trueGrandDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      } else if (dau2_truepdg == -211) {
        piMinus = daughter2;
        k0vtxpmtruepdg_val = dau2_truepdg;
        // Get process information for pi-
        if (daughter2->TrueObject) {
          AnaTrueParticlePD* trueDau2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
          if (trueDau2) {
            k0vtxpmtrueproc_val = static_cast<Int_t>(trueDau2->ProcessStart);
            k0vtxpmtrueendproc_val = static_cast<Int_t>(trueDau2->ProcessEnd);

            // Fill daughter PDGs for pi-
            for (int j = 0; j < trueDau2->Daughters.size() && j < 200; j++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t k = 0; k < nTrueParts; k++) {
                if (trueParticles[k] && trueParticles[k]->ID == trueDau2->Daughters[j]) {
                  AnaTrueParticlePD* trueGrandDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[k]);
                  if (trueGrandDaughter) {
                    k0vtxtruepmtruedaupdg_val[j] = trueGrandDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      // Calculate particle-specific variables for pi+
      if (piPlus) {
        k0vtxpplength_val = pdAnaUtils::ComputeTrackLengthFromHitPosition(piPlus);

        // Reconstructed information
        k0vtxpprecondau_val = piPlus->Daughters.size();
        k0vtxpprecostartdir_val[0] = piPlus->DirectionStart[0];
        k0vtxpprecostartdir_val[1] = piPlus->DirectionStart[1];
        k0vtxpprecostartdir_val[2] = piPlus->DirectionStart[2];
        k0vtxpprecoenddir_val[0] = piPlus->DirectionEnd[0];
        k0vtxpprecoenddir_val[1] = piPlus->DirectionEnd[1];
        k0vtxpprecoenddir_val[2] = piPlus->DirectionEnd[2];

        // Reco PDG of daughters
        for (int i = 0; i < piPlus->Daughters.size() && i < 200; i++) {
          AnaParticlePD* daughter = static_cast<AnaParticlePD*>(piPlus->Daughters[i]);
          if (daughter) {
            // Use plane 2 as default, or first non-default value
            if (daughter->ReconPDG[2] != -999) {
              k0vtxpprecopdgdau_val[i] = daughter->ReconPDG[2];
            } else if (daughter->ReconPDG[0] != -999) {
              k0vtxpprecopdgdau_val[i] = daughter->ReconPDG[0];
            } else if (daughter->ReconPDG[1] != -999) {
              k0vtxpprecopdgdau_val[i] = daughter->ReconPDG[1];
            }
          }
        }

        // True information
        if (piPlus->TrueObject) {
          AnaTrueParticlePD* truePiPlus = static_cast<AnaTrueParticlePD*>(piPlus->TrueObject);
          if (truePiPlus) {
            k0vtxpptruendau_val = truePiPlus->Daughters.size();
            k0vtxpptruestartdir_val[0] = truePiPlus->Direction[0];
            k0vtxpptruestartdir_val[1] = truePiPlus->Direction[1];
            k0vtxpptruestartdir_val[2] = truePiPlus->Direction[2];
            k0vtxpptrueenddir_val[0] = truePiPlus->DirectionEnd[0];
            k0vtxpptrueenddir_val[1] = truePiPlus->DirectionEnd[1];
            k0vtxpptrueenddir_val[2] = truePiPlus->DirectionEnd[2];

            // True PDG of daughters
            for (int i = 0; i < truePiPlus->Daughters.size() && i < 200; i++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t j = 0; j < nTrueParts; j++) {
                if (trueParticles[j] && trueParticles[j]->ID == truePiPlus->Daughters[i]) {
                  AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[j]);
                  if (trueDaughter) {
                    k0vtxpptruepdgdau_val[i] = trueDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      // Calculate particle-specific variables for pi-
      if (piMinus) {
        k0vtxpmlength_val = pdAnaUtils::ComputeTrackLengthFromHitPosition(piMinus);

        // Reconstructed information
        k0vtxpmrecondau_val = piMinus->Daughters.size();
        k0vtxpmrecostartdir_val[0] = piMinus->DirectionStart[0];
        k0vtxpmrecostartdir_val[1] = piMinus->DirectionStart[1];
        k0vtxpmrecostartdir_val[2] = piMinus->DirectionStart[2];
        k0vtxpmrecoenddir_val[0] = piMinus->DirectionEnd[0];
        k0vtxpmrecoenddir_val[1] = piMinus->DirectionEnd[1];
        k0vtxpmrecoenddir_val[2] = piMinus->DirectionEnd[2];

        // Reco PDG of daughters
        for (int i = 0; i < piMinus->Daughters.size() && i < 200; i++) {
          // std::cout << "Daughter " << i << " with ID " << piMinus->Daughters[i] << std::endl;
          AnaParticlePD* daughter = static_cast<AnaParticlePD*>(piMinus->Daughters[i]);
          if (daughter) {
            // Use plane 2 as default, or first non-default value
            if (daughter->ReconPDG[2] != -999) {
              k0vtxpmrecopdgdau_val[i] = daughter->ReconPDG[2];
            } else if (daughter->ReconPDG[0] != -999) {
              k0vtxpmrecopdgdau_val[i] = daughter->ReconPDG[0];
            } else if (daughter->ReconPDG[1] != -999) {
              k0vtxpmrecopdgdau_val[i] = daughter->ReconPDG[1];
            }
            // std::cout << "Daughter " << i << " reco PDG: " << k0vtxpmrecopdgdau_val[i] << std::endl;
          }
        }

        // True information
        if (piMinus->TrueObject) {
          AnaTrueParticlePD* truePiMinus = static_cast<AnaTrueParticlePD*>(piMinus->TrueObject);
          if (truePiMinus) {
            k0vtxpmtruendau_val = truePiMinus->Daughters.size();
            k0vtxpmtruestartdir_val[0] = truePiMinus->Direction[0];
            k0vtxpmtruestartdir_val[1] = truePiMinus->Direction[1];
            k0vtxpmtruestartdir_val[2] = truePiMinus->Direction[2];
            k0vtxpmtrueenddir_val[0] = truePiMinus->DirectionEnd[0];
            k0vtxpmtrueenddir_val[1] = truePiMinus->DirectionEnd[1];
            k0vtxpmtrueenddir_val[2] = truePiMinus->DirectionEnd[2];

            // True PDG of daughters
            for (int i = 0; i < truePiMinus->Daughters.size() && i < 200; i++) {
              // Find the true particle corresponding to this daughter ID
              AnaTrueParticleB** trueParticles = event.TrueParticles;
              Int_t nTrueParts = event.nTrueParticles;
              for (Int_t j = 0; j < nTrueParts; j++) {
                if (trueParticles[j] && trueParticles[j]->ID == truePiMinus->Daughters[i]) {
                  // std::cout << "True particle found for daughter " << i << " with ID " << truePiMinus->Daughters[i] << std::endl;
                  AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(trueParticles[j]);
                  if (trueDaughter) {
                    // std::cout << "True daughter PDG: " << trueDaughter->PDG << std::endl;
                    k0vtxpmtruepdgdau_val[i] = trueDaughter->PDG;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      // ========================================
      // RECONSTRUCTED VERTEX SYSTEM VARIABLES
      // ========================================

      // Calculate reconstructed momentum vectors
      Float_t mom1_reco[3] = {-999.0, -999.0, -999.0};
      Float_t mom2_reco[3] = {-999.0, -999.0, -999.0};
      Float_t mom1_mag = -999.0;
      Float_t mom2_mag = -999.0;

      // Get kinetic energies
      Float_t kinetic_energy1 = pdAnaUtils::ComputeKineticEnergy(*daughter1);
      Float_t kinetic_energy2 = pdAnaUtils::ComputeKineticEnergy(*daughter2);
      Float_t deposited_energy1 = pdAnaUtils::ComputeDepositedEnergy(daughter1);
      Float_t deposited_energy2 = pdAnaUtils::ComputeDepositedEnergy(daughter2);

      // Calculate momentum using the same method as mass calculation
      const Float_t pion_mass = 139.57; // MeV/c²
      if (kinetic_energy1 > 0 && kinetic_energy2 > 0) {
        // Method 1: Use kinetic energy
        Float_t KE1_MeV = kinetic_energy1 * 1000.0;
        Float_t KE2_MeV = kinetic_energy2 * 1000.0;
        Float_t E1 = KE1_MeV + pion_mass;
        Float_t E2 = KE2_MeV + pion_mass;
        mom1_mag = sqrt(E1*E1 - pion_mass*pion_mass);
        mom2_mag = sqrt(E2*E2 - pion_mass*pion_mass);
      } else if (deposited_energy1 > 0 && deposited_energy2 > 0) {
        // Method 2: Use deposited energy
        Float_t E1 = deposited_energy1 + pion_mass;
        Float_t E2 = deposited_energy2 + pion_mass;
        mom1_mag = sqrt(E1*E1 - pion_mass*pion_mass);
        mom2_mag = sqrt(E2*E2 - pion_mass*pion_mass);
      } else {
        // Method 3: Use range-momentum relation
        Float_t p1_mag = pdAnaUtils::ComputeRangeMomentum(daughter1->Length, 211);
        Float_t p2_mag = pdAnaUtils::ComputeRangeMomentum(daughter2->Length, 211);
        mom1_mag = p1_mag * 1000.0; // Convert GeV/c to MeV/c
        mom2_mag = p2_mag * 1000.0; // Convert GeV/c to MeV/c
      }

      // Calculate momentum vectors
      mom1_reco[0] = mom1_mag * daughter1->DirectionStart[0];
      mom1_reco[1] = mom1_mag * daughter1->DirectionStart[1];
      mom1_reco[2] = mom1_mag * daughter1->DirectionStart[2];
      mom2_reco[0] = mom2_mag * daughter2->DirectionStart[0];
      mom2_reco[1] = mom2_mag * daughter2->DirectionStart[1];
      mom2_reco[2] = mom2_mag * daughter2->DirectionStart[2];

      // Calculate total system momentum
      Float_t total_mom_reco[3] = {mom1_reco[0] + mom2_reco[0],
                                   mom1_reco[1] + mom2_reco[1],
                                   mom1_reco[2] + mom2_reco[2]};
      k0vtxrecomom_val = sqrt(total_mom_reco[0]*total_mom_reco[0] +
                              total_mom_reco[1]*total_mom_reco[1] +
                              total_mom_reco[2]*total_mom_reco[2]);

      // Calculate system direction (unit vector)
      if (k0vtxrecomom_val > 0) {
        k0vtxrecodir_val[0] = total_mom_reco[0] / k0vtxrecomom_val;
        k0vtxrecodir_val[1] = total_mom_reco[1] / k0vtxrecomom_val;
        k0vtxrecodir_val[2] = total_mom_reco[2] / k0vtxrecomom_val;
      }

      // Calculate opening angle between daughters
      if (mom1_mag > 0 && mom2_mag > 0) {
        Float_t dot_product = mom1_reco[0]*mom2_reco[0] + mom1_reco[1]*mom2_reco[1] + mom1_reco[2]*mom2_reco[2];
        Float_t cos_angle = dot_product / (mom1_mag * mom2_mag);
        if (cos_angle >= -1.0 && cos_angle <= 1.0) {
          k0vtxrecoopening_val = acos(cos_angle) * 180.0 / M_PI; // Convert to degrees
        }
      }

      // Calculate system energy
      Float_t E1_reco = sqrt(mom1_mag*mom1_mag + pion_mass*pion_mass);
      Float_t E2_reco = sqrt(mom2_mag*mom2_mag + pion_mass*pion_mass);
      k0vtxrecoenergy_val = E1_reco + E2_reco;

      // Calculate system mass (invariant mass)
      Float_t total_energy_reco = k0vtxrecoenergy_val;
      Float_t total_mom_sq_reco = k0vtxrecomom_val * k0vtxrecomom_val;
      Float_t mass_sq_reco = total_energy_reco * total_energy_reco - total_mom_sq_reco;
      if (mass_sq_reco >= 0) {
        k0vtxrecomass_val = sqrt(mass_sq_reco);
        k0vtxrecoinvariantmass_val = k0vtxrecomass_val;
      }

      // Calculate ctau (if we have vertex position and system momentum)
      if (candidate->Parent) {
        Float_t dx = candidate->PositionStart[0] - candidate->Parent->PositionEnd[0];
        Float_t dy = candidate->PositionStart[1] - candidate->Parent->PositionEnd[1];
        Float_t dz = candidate->PositionStart[2] - candidate->Parent->PositionEnd[2];
        Float_t distance = sqrt(dx*dx + dy*dy + dz*dz);
        if (k0vtxrecomom_val > 0) {
          k0vtxrecoctau_val = distance * pion_mass / k0vtxrecomom_val; // ctau = distance * m/p
        }
      }

      // Calculate angle with beam direction (if available)
      AnaBeamPD* beam = static_cast<AnaBeamPD*>(eventPD->Beam);
      if (beam && beam->BeamParticle) {
        AnaParticlePD* beamParticle = static_cast<AnaParticlePD*>(beam->BeamParticle);
        if (beamParticle && beamParticle->DirectionStart[0] != -999) {
          Float_t beam_dot_system = beamParticle->DirectionStart[0] * k0vtxrecodir_val[0] +
                                    beamParticle->DirectionStart[1] * k0vtxrecodir_val[1] +
                                    beamParticle->DirectionStart[2] * k0vtxrecodir_val[2];
          if (beam_dot_system >= -1.0 && beam_dot_system <= 1.0) {
            k0vtxrecoangle_val = acos(beam_dot_system) * 180.0 / M_PI; // Convert to degrees
          }
        }
      }

      // ========================================
      // COUNT POTENTIAL PARENT PARTICLES
      // ========================================

      // Count particles whose end positions are within DaughterDistance of either vertex particle's start position
      // Exclude: the parent particle and the two vertex particles themselves
      const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");

      // Get start positions of both vertex particles
      TVector3 vtxPart1Pos(daughter1->PositionStart[0], daughter1->PositionStart[1], daughter1->PositionStart[2]);
      TVector3 vtxPart2Pos(daughter2->PositionStart[0], daughter2->PositionStart[1], daughter2->PositionStart[2]);

      // Loop over all reconstructed particles in the event
      for (int i = 0; i < eventPD->nParticles; i++) {
        AnaParticlePD* particle = static_cast<AnaParticlePD*>(eventPD->Particles[i]);
        if (!particle) continue;

        // Skip the parent particle
        if (candidate->Parent && particle == candidate->Parent) continue;

        // Skip the two vertex particles
        if (particle == daughter1 || particle == daughter2) continue;

        // Check if particle's end position is within DaughterDistance of either vertex particle
        TVector3 particleEndPos(particle->PositionEnd[0], particle->PositionEnd[1], particle->PositionEnd[2]);

        // Skip if end position is invalid
        if (particle->PositionEnd[0] < -900) continue;

        // Calculate distance to both vertex particles
        double dist1 = (particleEndPos - vtxPart1Pos).Mag();
        double dist2 = (particleEndPos - vtxPart2Pos).Mag();

        // Count if within threshold of either vertex particle
        if (dist1 <= maxDaughterDistance || dist2 <= maxDaughterDistance) {
          k0vtxnpotpar_val++;
        }
      }

      // ========================================
      // TRUE VERTEX SYSTEM VARIABLES
      // ========================================

      if (daughter1->TrueObject && daughter2->TrueObject) {
      AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
      AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

      if (trueDaughter1 && trueDaughter2) {
        // Calculate true momentum vectors
        Float_t mom1_true[3] = {trueDaughter1->Momentum * trueDaughter1->Direction[0],
                                trueDaughter1->Momentum * trueDaughter1->Direction[1],
                                trueDaughter1->Momentum * trueDaughter1->Direction[2]};
        Float_t mom2_true[3] = {trueDaughter2->Momentum * trueDaughter2->Direction[0],
                                trueDaughter2->Momentum * trueDaughter2->Direction[1],
                                trueDaughter2->Momentum * trueDaughter2->Direction[2]};

        // Calculate total system momentum
        Float_t total_mom_true[3] = {mom1_true[0] + mom2_true[0],
                                     mom1_true[1] + mom2_true[1],
                                     mom1_true[2] + mom2_true[2]};
        k0vtxtruemom_val = sqrt(total_mom_true[0]*total_mom_true[0] +
                                total_mom_true[1]*total_mom_true[1] +
                                total_mom_true[2]*total_mom_true[2]);

        // Calculate system direction (unit vector)
        if (k0vtxtruemom_val > 0) {
          k0vtxtruedir_val[0] = total_mom_true[0] / k0vtxtruemom_val;
          k0vtxtruedir_val[1] = total_mom_true[1] / k0vtxtruemom_val;
          k0vtxtruedir_val[2] = total_mom_true[2] / k0vtxtruemom_val;
        }

        // Calculate opening angle between daughters
        if (trueDaughter1->Momentum > 0 && trueDaughter2->Momentum > 0) {
          Float_t dot_product = mom1_true[0]*mom2_true[0] + mom1_true[1]*mom2_true[1] + mom1_true[2]*mom2_true[2];
          Float_t cos_angle = dot_product / (trueDaughter1->Momentum * trueDaughter2->Momentum);
          if (cos_angle >= -1.0 && cos_angle <= 1.0) {
            k0vtxtrueopening_val = acos(cos_angle) * 180.0 / M_PI; // Convert to degrees
          }
        }

        // Calculate system energy
        Float_t mass1_true = anaUtils::GetParticleMass(ParticleId::GetParticle(trueDaughter1->PDG));
        Float_t mass2_true = anaUtils::GetParticleMass(ParticleId::GetParticle(trueDaughter2->PDG));
        Float_t E1_true = sqrt(trueDaughter1->Momentum*trueDaughter1->Momentum + mass1_true*mass1_true);
        Float_t E2_true = sqrt(trueDaughter2->Momentum*trueDaughter2->Momentum + mass2_true*mass2_true);
        k0vtxtrueenergy_val = E1_true + E2_true;

        // Calculate system mass (invariant mass)
        Float_t total_energy_true = k0vtxtrueenergy_val;
        Float_t total_mom_sq_true = k0vtxtruemom_val * k0vtxtruemom_val;
        Float_t mass_sq_true = total_energy_true * total_energy_true - total_mom_sq_true;
        if (mass_sq_true >= 0) {
          k0vtxtruemass_val = sqrt(mass_sq_true);
          k0vtxtrueinvariantmass_val = k0vtxtruemass_val;
        }

        // Calculate ctau (if we have true positions)
        Float_t dx_true = trueDaughter1->Position[0] - candidate->Parent->PositionEnd[0];
        Float_t dy_true = trueDaughter1->Position[1] - candidate->Parent->PositionEnd[1];
        Float_t dz_true = trueDaughter1->Position[2] - candidate->Parent->PositionEnd[2];
        Float_t distance_true = sqrt(dx_true*dx_true + dy_true*dy_true + dz_true*dz_true);
        if (k0vtxtruemom_val > 0) {
          k0vtxtruectau_val = distance_true * mass1_true / k0vtxtruemom_val; // ctau = distance * m/p
        }

        // Calculate angle with beam direction (if available)
        AnaBeamPD* beam = static_cast<AnaBeamPD*>(eventPD->Beam);
        if (beam && beam->BeamParticle && beam->BeamParticle->TrueObject) {
          AnaTrueParticlePD* beamTrueParticle = static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);
          if (beamTrueParticle) {
            Float_t beam_dot_system = beamTrueParticle->Direction[0] * k0vtxtruedir_val[0] +
                                      beamTrueParticle->Direction[1] * k0vtxtruedir_val[1] +
                                      beamTrueParticle->Direction[2] * k0vtxtruedir_val[2];
            if (beam_dot_system >= -1.0 && beam_dot_system <= 1.0) {
              k0vtxtrueangle_val = acos(beam_dot_system) * 180.0 / M_PI; // Convert to degrees
            }
          }
        }
      }
      }
    }
  }

  // Fill vertex system variables
  output.FillVectorVar(k0vtxrecomom, k0vtxrecomom_val);
  output.FillVectorVar(k0vtxtruemom, k0vtxtruemom_val);
  output.FillVectorVar(k0vtxrecoenergy, k0vtxrecoenergy_val);
  output.FillVectorVar(k0vtxtrueenergy, k0vtxtrueenergy_val);
  output.FillVectorVar(k0vtxrecomass, k0vtxrecomass_val);
  output.FillVectorVar(k0vtxtruemass, k0vtxtruemass_val);
  output.FillMatrixVarFromArray(k0vtxrecodir, k0vtxrecodir_val, 3);
  output.FillMatrixVarFromArray(k0vtxtruedir, k0vtxtruedir_val, 3);
  output.FillVectorVar(k0vtxrecoopening, k0vtxrecoopening_val);
  output.FillVectorVar(k0vtxtrueopening, k0vtxtrueopening_val);
  output.FillVectorVar(k0vtxrecoinvariantmass, k0vtxrecoinvariantmass_val);
  output.FillVectorVar(k0vtxtrueinvariantmass, k0vtxtrueinvariantmass_val);
  output.FillVectorVar(k0vtxrecoctau, k0vtxrecoctau_val);
  output.FillVectorVar(k0vtxtruectau, k0vtxtruectau_val);
  output.FillVectorVar(k0vtxrecoangle, k0vtxrecoangle_val);
  output.FillVectorVar(k0vtxtrueangle, k0vtxtrueangle_val);
  output.FillVectorVar(k0vtxnpotpar, k0vtxnpotpar_val);
  output.FillVectorVar(k0vtxpplength, k0vtxpplength_val);
  output.FillVectorVar(k0vtxpmlength, k0vtxpmlength_val);
  output.FillVectorVar(k0vtxpptruepdg, k0vtxpptruepdg_val);
  output.FillVectorVar(k0vtxpmtruepdg, k0vtxpmtruepdg_val);
  output.FillVectorVar(k0vtxpptruendau, k0vtxpptruendau_val);
  output.FillVectorVar(k0vtxpmtruendau, k0vtxpmtruendau_val);
  output.FillVectorVar(k0vtxpprecondau, k0vtxpprecondau_val);
  output.FillVectorVar(k0vtxpmrecondau, k0vtxpmrecondau_val);

  // Fill daughter PDG arrays (using matrix variables)
  output.FillMatrixVarFromArray(k0vtxpptruepdgdau, k0vtxpptruepdgdau_val, 200);
  output.FillMatrixVarFromArray(k0vtxpmtruepdgdau, k0vtxpmtruepdgdau_val, 200);
  output.FillMatrixVarFromArray(k0vtxpprecopdgdau, k0vtxpprecopdgdau_val, 200);
  output.FillMatrixVarFromArray(k0vtxpmrecopdgdau, k0vtxpmrecopdgdau_val, 200);

  // Fill direction arrays (3D matrix variables)
  output.FillMatrixVarFromArray(k0vtxpptruestartdir, k0vtxpptruestartdir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpmtruestartdir, k0vtxpmtruestartdir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpptrueenddir, k0vtxpptrueenddir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpmtrueenddir, k0vtxpmtrueenddir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpprecostartdir, k0vtxpprecostartdir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpmrecostartdir, k0vtxpmrecostartdir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpprecoenddir, k0vtxpprecoenddir_val, 3);
  output.FillMatrixVarFromArray(k0vtxpmrecoenddir, k0vtxpmrecoenddir_val, 3);
  output.FillVectorVar(k0vtxpptrueproc, k0vtxpptrueproc_val);
  output.FillVectorVar(k0vtxpmtrueproc, k0vtxpmtrueproc_val);
  output.FillVectorVar(k0vtxpptrueendproc, k0vtxpptrueendproc_val);
  output.FillVectorVar(k0vtxpmtrueendproc, k0vtxpmtrueendproc_val);
  output.FillMatrixVarFromArray(k0vtxtruepptruedaupdg, k0vtxtruepptruedaupdg_val, 200);
  output.FillMatrixVarFromArray(k0vtxtruepmtruedaupdg, k0vtxtruepmtruedaupdg_val, 200);

  // ========================================
  // PROTON BROTHERS OF K0 (daughters of K0 parent)
  // ========================================
  // Loop over reconstructed daughters of the K0 parent and identify protons
  if (candidate->Parent) {
    for (int i = 0; i < candidate->Parent->Daughters.size(); i++) {
      AnaParticlePD* daughter = static_cast<AnaParticlePD*>(candidate->Parent->Daughters[i]);
      if (daughter && daughter->TrueObject) {
        AnaTrueParticlePD* trueDaughter = static_cast<AnaTrueParticlePD*>(daughter->TrueObject);
        if (trueDaughter && trueDaughter->PDG == 2212) { // Proton PDG code
          k0brothnproton_val++;

          // Fill variables for the first proton found
          if (k0brothnproton_val == 1) {
            // True information
            k0brothtruepdg_val = trueDaughter->PDG;
            k0brothtrueproc_val = static_cast<Int_t>(trueDaughter->ProcessStart);
            k0brothtrueendproc_val = static_cast<Int_t>(trueDaughter->ProcessEnd);
            k0brothtruendau_val = trueDaughter->Daughters.size();

            // Reconstructed information
            k0brothrecondau_val = daughter->Daughters.size();
            k0brothlength_val = pdAnaUtils::ComputeTrackLengthFromHitPosition(daughter);

            // Reconstructed PDG - use plane 2 as default, or first non-default value
            if (daughter->ReconPDG[2] != -999) {
              k0brothrecopdg_val = daughter->ReconPDG[2];
            } else if (daughter->ReconPDG[0] != -999) {
              k0brothrecopdg_val = daughter->ReconPDG[0];
            } else if (daughter->ReconPDG[1] != -999) {
              k0brothrecopdg_val = daughter->ReconPDG[1];
            }
          }
        }
      }
    }
  }

  // Fill brother variables
  output.FillVectorVar(k0brothnproton, k0brothnproton_val);
  output.FillVectorVar(k0brothtruepdg, k0brothtruepdg_val);
  output.FillVectorVar(k0brothrecopdg, k0brothrecopdg_val);
  output.FillVectorVar(k0brothtrueproc, k0brothtrueproc_val);
  output.FillVectorVar(k0brothtrueendproc, k0brothtrueendproc_val);
  output.FillVectorVar(k0brothtruendau, k0brothtruendau_val);
  output.FillVectorVar(k0brothrecondau, k0brothrecondau_val);
  output.FillVectorVar(k0brothlength, k0brothlength_val);

  // ========================================
  // 4. FILL K0 DAUGHTER VARIABLES (loop already implemented)
  // ========================================

  // Note: Daughter variables are filled in the FillNeutralKaonVariables_Candidates function
  // which loops through all candidates and fills daughter variables for each one

  // output.IncrementCounter(nk0);
}
