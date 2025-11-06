#include "neutralKaonEventDisplay.hxx"
#include "ToyBoxNeutralKaon.hxx"
#include "pdDataClasses.hxx"
#include "OutputManager.hxx"
#include <TEveGeoShape.h>
#include <TGeoSphere.h>
#include <TGeoVolume.h>
#include <TEvePointSet.h>
#include <TEveLine.h>
#include <TEveScene.h>
#include <TEveTrans.h>
#include <TTree.h>
#include <iostream>

//********************************************************************
neutralKaonEventDisplay::neutralKaonEventDisplay() : pdEventDisplay() {
//********************************************************************
    // Set class name for auto-detection
    _eventDisplayClassName = "neutralKaonEventDisplay";

    _nK0Candidates = 0;

    // Initialize arrays
    for (Int_t i = 0; i < kMaxK0; i++) {
        _k0_daughter1ID[i] = -1;
        _k0_daughter2ID[i] = -1;
        _k0_parentID[i] = -1;
        _k0_creationVtxRadius[i] = 0;
        _k0_annihilationVtxRadius[i] = 0;
        for (Int_t j = 0; j < 3; j++) {
            _k0_creationVtxPos[i][j] = -999;
            _k0_annihilationVtxPos[i][j] = -999;
            _k0_startPos[i][j] = -999;
            _k0_endPos[i][j] = -999;
        }
    }
}

//********************************************************************
neutralKaonEventDisplay::~neutralKaonEventDisplay() {
//********************************************************************
}

//********************************************************************
void neutralKaonEventDisplay::AddAnalysisVariables(OutputManager& output, Int_t tree_index) {
//********************************************************************
    // Add K0 candidate variables
    output.AddVectorVar(tree_index, edk0_daughter1ID, "ED_k0_daughter1ID", "I", "K0 daughter 1 unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_daughter2ID, "ED_k0_daughter2ID", "I", "K0 daughter 2 unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_parentID, "ED_k0_parentID", "I", "K0 parent unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_creationVtxRadius, "ED_k0_creationVtxRadius", "F", "K0 creation vertex radius", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_annihilationVtxRadius, "ED_k0_annihilationVtxRadius", "F", "K0 annihilation vertex radius", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_creationVtxPos, "ED_k0_creationVtxPos", "F", "K0 creation vertex positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annihilationVtxPos, "ED_k0_annihilationVtxPos", "F", "K0 annihilation vertex positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_startPos, "ED_k0_startPos", "F", "K0 start positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_endPos, "ED_k0_endPos", "F", "K0 end positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);

    std::cout << "neutralKaonEventDisplay::AddAnalysisVariables() - K0 variables added" << std::endl;
}

//********************************************************************
void neutralKaonEventDisplay::FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
    // Cast box to ToyBoxNeutralKaon
    const ToyBoxNeutralKaon* k0Box = static_cast<const ToyBoxNeutralKaon*>(box);

    // Fill K0 candidate data from neutral particle candidates
    if (k0Box && k0Box->neutralParticleCandidates.size() > 0) {
        for (size_t i = 0; i < k0Box->neutralParticleCandidates.size() && i < (size_t)kMaxK0; i++) {
            AnaNeutralParticlePD* neutralParticle = k0Box->neutralParticleCandidates[i];
            if (!neutralParticle) continue;

            // Get daughter particles from annihilation vertex
            Int_t daughter1ID = -1;
            Int_t daughter2ID = -1;
            if (neutralParticle->AnnihilationVertex && neutralParticle->AnnihilationVertex->Particles.size() >= 2) {
                daughter1ID = neutralParticle->AnnihilationVertex->Particles[0]->UniqueID;
                daughter2ID = neutralParticle->AnnihilationVertex->Particles[1]->UniqueID;
            }

            // Get parent ID
            Int_t parentID = neutralParticle->Parent ? neutralParticle->Parent->UniqueID : -1;

            // Fill K0 data
            output.FillVectorVar(edk0_daughter1ID, daughter1ID);
            output.FillVectorVar(edk0_daughter2ID, daughter2ID);
            output.FillVectorVar(edk0_parentID, parentID);

            // Vertex radii (use a default if not available)
            Float_t creationRadius = 30.0; // Default radius
            Float_t annihilationRadius = 30.0;

            output.FillVectorVar(edk0_creationVtxRadius, creationRadius);
            output.FillVectorVar(edk0_annihilationVtxRadius, annihilationRadius);

            // Creation vertex position (from parent end position or neutral particle position)
            Float_t creationPos[3] = {-999, -999, -999};
            if (neutralParticle->CreationVertex) {
                creationPos[0] = (Float_t)neutralParticle->CreationVertex->Position[0];
                creationPos[1] = (Float_t)neutralParticle->CreationVertex->Position[1];
                creationPos[2] = (Float_t)neutralParticle->CreationVertex->Position[2];
            } else if (neutralParticle->Parent) {
                creationPos[0] = (Float_t)neutralParticle->Parent->PositionEnd[0];
                creationPos[1] = (Float_t)neutralParticle->Parent->PositionEnd[1];
                creationPos[2] = (Float_t)neutralParticle->Parent->PositionEnd[2];
            }

            // Annihilation vertex position
            Float_t annihilationPos[3] = {-999, -999, -999};
            if (neutralParticle->AnnihilationVertex) {
                annihilationPos[0] = (Float_t)neutralParticle->AnnihilationVertex->Position[0];
                annihilationPos[1] = (Float_t)neutralParticle->AnnihilationVertex->Position[1];
                annihilationPos[2] = (Float_t)neutralParticle->AnnihilationVertex->Position[2];
            }

            // K0 trajectory (from creation to annihilation)
            Float_t startPos[3] = {creationPos[0], creationPos[1], creationPos[2]};
            Float_t endPos[3] = {annihilationPos[0], annihilationPos[1], annihilationPos[2]};

            output.FillMatrixVarFromArray(edk0_creationVtxPos, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annihilationVtxPos, annihilationPos, 3);
            output.FillMatrixVarFromArray(edk0_startPos, startPos, 3);
            output.FillMatrixVarFromArray(edk0_endPos, endPos, 3);

            // Increment K0 counter
            output.IncrementCounter(ednK0Candidates);
        }
    }
}

//********************************************************************
bool neutralKaonEventDisplay::ReadAnalysisData(TTree* tree) {
//********************************************************************
    // Set branch addresses for K0 candidate data
    tree->SetBranchAddress("ED_nK0Candidates", &_nK0Candidates);
    tree->SetBranchAddress("ED_k0_daughter1ID", _k0_daughter1ID);
    tree->SetBranchAddress("ED_k0_daughter2ID", _k0_daughter2ID);
    tree->SetBranchAddress("ED_k0_parentID", _k0_parentID);
    tree->SetBranchAddress("ED_k0_creationVtxRadius", _k0_creationVtxRadius);
    tree->SetBranchAddress("ED_k0_annihilationVtxRadius", _k0_annihilationVtxRadius);
    tree->SetBranchAddress("ED_k0_creationVtxPos", _k0_creationVtxPos);
    tree->SetBranchAddress("ED_k0_annihilationVtxPos", _k0_annihilationVtxPos);
    tree->SetBranchAddress("ED_k0_startPos", _k0_startPos);
    tree->SetBranchAddress("ED_k0_endPos", _k0_endPos);

    // Read the current entry
    tree->GetEntry(tree->GetReadEntry());

    std::cout << "Read K0 analysis data: " << _nK0Candidates << " K0 candidates" << std::endl;

    return true;
}

//********************************************************************
void neutralKaonEventDisplay::DrawAnalysisContent3D(TEveScene* scene) {
//********************************************************************
    // Draw K0 candidates
    for (Int_t i = 0; i < _nK0Candidates && i < kMaxK0; i++) {
        Int_t parentColor = kBlue;

        // Creation vertex sphere
        Float_t creationX = _k0_creationVtxPos[i][0];
        Float_t creationY = _k0_creationVtxPos[i][1];
        Float_t creationZ = _k0_creationVtxPos[i][2];

        if (creationX > -900 && creationY > -900 && creationZ > -900) {
            Float_t radius = _k0_creationVtxRadius[i];
            TGeoSphere* creationSphere = new TGeoSphere(0, radius);
            TGeoVolume* creationVol = new TGeoVolume(Form("CreationSphere_%d", i), creationSphere);
            TEveGeoShape* creationShape = new TEveGeoShape(Form("K0 #%d Creation Radius (%.0f cm)", i, radius));
            creationShape->SetShape(creationVol->GetShape());
            creationShape->SetMainColor(parentColor);
            creationShape->SetMainTransparency(80);
            TEveTrans& creationTrans = creationShape->RefMainTrans();
            creationTrans.SetPos(creationX, creationY, creationZ);
            scene->AddElement(creationShape);

            // Creation vertex marker
            TEvePointSet* vtx = new TEvePointSet(Form("K0 #%d Creation Vertex", i));
            vtx->SetNextPoint(creationX, creationY, creationZ);
            vtx->SetMarkerStyle(29); // Star
            vtx->SetMarkerSize(3.0);
            vtx->SetMainColor(parentColor);
            scene->AddElement(vtx);
        }

        // Annihilation vertex sphere
        Float_t annihilationX = _k0_annihilationVtxPos[i][0];
        Float_t annihilationY = _k0_annihilationVtxPos[i][1];
        Float_t annihilationZ = _k0_annihilationVtxPos[i][2];

        if (annihilationX > -900 && annihilationY > -900 && annihilationZ > -900) {
            Float_t annihilationRadius = _k0_annihilationVtxRadius[i];
            TGeoSphere* annihilationSphere = new TGeoSphere(0, annihilationRadius);
            TGeoVolume* annihilationVol = new TGeoVolume(Form("AnnihilationSphere_%d", i), annihilationSphere);
            TEveGeoShape* annihilationShape = new TEveGeoShape(Form("K0 #%d Annihilation Radius (%.0f cm)", i, annihilationRadius));
            annihilationShape->SetShape(annihilationVol->GetShape());
            annihilationShape->SetMainColor(kRed);
            annihilationShape->SetMainTransparency(80);
            TEveTrans& annihilationTrans = annihilationShape->RefMainTrans();
            annihilationTrans.SetPos(annihilationX, annihilationY, annihilationZ);
            scene->AddElement(annihilationShape);

            // Annihilation vertex marker
            TEvePointSet* vtx = new TEvePointSet(Form("K0 #%d Decay Vertex", i));
            vtx->SetNextPoint(annihilationX, annihilationY, annihilationZ);
            vtx->SetMarkerStyle(29); // Star
            vtx->SetMarkerSize(3.0);
            vtx->SetMainColor(kRed);
            scene->AddElement(vtx);
        }

        // K0 trajectory (neutral particle path - magenta dashed)
        if (_k0_startPos[i][0] > -900 && _k0_endPos[i][0] > -900) {
            TEveLine* k0 = new TEveLine(Form("K0 #%d Trajectory", i));
            k0->SetPoint(0, _k0_startPos[i][0], _k0_startPos[i][1], _k0_startPos[i][2]);
            k0->SetPoint(1, _k0_endPos[i][0], _k0_endPos[i][1], _k0_endPos[i][2]);
            k0->SetMainColor(kMagenta);
            k0->SetLineWidth(3);
            k0->SetLineStyle(2); // Dashed
            scene->AddElement(k0);
        }
    }

    if (_nK0Candidates > 0) {
        std::cout << "Drew " << _nK0Candidates << " K0 candidates in 3D" << std::endl;
    }
}

//********************************************************************
void neutralKaonEventDisplay::DrawAnalysisContent2D(TEveProjectionManager* manager, const std::string& projection_type) {
//********************************************************************
    // K0 vertices will be automatically projected from 3D via ImportElements
    // No additional 2D-specific drawing needed
    (void)manager;
    (void)projection_type;
}

