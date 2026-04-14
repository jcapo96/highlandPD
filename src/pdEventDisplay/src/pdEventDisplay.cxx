#include "pdEventDisplay.hxx"
#include "OutputManager.hxx"
#include "BasicUtils.hxx"
#include "ToyBoxNeutralKaon.hxx"
#include <TEveManager.h>
#include <TEveLine.h>
#include <TEvePointSet.h>
#include <TEveScene.h>
#include <TEveGeoShape.h>
#include <TEveElement.h>
#include <TEveTrans.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoSphere.h>
#include <TEveProjectionManager.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>
#include <set>
#include <map>

namespace {
Float_t HitMarkerSizeFromdEdx(Float_t dEdx) {
    // Scale hit marker size with local dEdx to enhance stopping-particle visibility.
    if (!std::isfinite(dEdx) || dEdx < 0.f) return 0.45f;
    const Float_t capped = std::max(0.f, std::min(10.f, dEdx));
    return 0.45f + 0.20f * capped; // [0.45, 2.45]
}

Int_t HitSizeBinFromdEdx(Float_t dEdx) {
    const Float_t size = HitMarkerSizeFromdEdx(dEdx);
    if (size < 0.8f) return 0;
    if (size < 1.2f) return 1;
    if (size < 1.8f) return 2;
    return 3;
}

Float_t HitSizeForBin(Int_t bin) {
    switch (bin) {
        case 0: return 0.55f;
        case 1: return 0.95f;
        case 2: return 1.45f;
        default: return 2.15f;
    }
}
}

//********************************************************************
pdEventDisplay::pdEventDisplay() : EventDisplayBase() {
//********************************************************************
    // Set class name for auto-detection
    _eventDisplayClassName = "pdEventDisplay";

    InitializeParticleColors();

    // Allocate hit arrays
    _hit_X = new Float_t[kMaxHits];
    _hit_Y = new Float_t[kMaxHits];
    _hit_Z = new Float_t[kMaxHits];
    _hit_dEdx = new Float_t[kMaxHits];

    // Initialize data members
    _nParticles = 0;
    _totalHits = 0;
    _nK0Candidates = 0;
    _hasTrueK0 = false;
}

//********************************************************************
pdEventDisplay::~pdEventDisplay() {
//********************************************************************
    delete[] _hit_X;
    delete[] _hit_Y;
    delete[] _hit_Z;
    delete[] _hit_dEdx;
}

//********************************************************************
void pdEventDisplay::InitializeParticleColors() {
//********************************************************************
    // Initialize particle colors based on PDG codes
    _particleColors[2212] = kOrange;     // Proton
    _particleColors[-2212] = kOrange;    // Anti-proton
    _particleColors[211] = kBlue;        // Pi+
    _particleColors[-211] = kRed;        // Pi-
    _particleColors[111] = kCyan;        // Pi0
    _particleColors[321] = kGreen;       // K+
    _particleColors[-321] = kYellow;     // K-
    _particleColors[310] = kMagenta;     // K0_S
    _particleColors[130] = kMagenta;     // K0_L
    _particleColors[311] = kMagenta;     // K0
    _particleColors[13] = kYellow;       // Mu-
    _particleColors[-13] = kYellow;      // Mu+
    _particleColors[11] = kTeal;         // e-
    _particleColors[-11] = kTeal;        // e+
    _particleColors[22] = kGray+2;       // gamma
    _particleColors[2112] = kAzure+1;    // Neutron
    _particleColors[-2112] = kAzure+1;   // Anti-neutron
}

//********************************************************************
int pdEventDisplay::GetParticleColor(int pdg) {
//********************************************************************
    if (_particleColors.find(pdg) != _particleColors.end()) {
        return _particleColors[pdg];
    }
    return kBlack; // Default color for unknown particles
}

//********************************************************************
std::string pdEventDisplay::GetParticleName(int pdg) {
//********************************************************************
    switch(abs(pdg)) {
        case 2212: return "Proton";
        case 211: return (pdg > 0) ? "π+" : "π-";
        case 111: return "π0";
        case 321: return (pdg > 0) ? "K+" : "K-";
        case 310: return "K0_S";
        case 130: return "K0_L";
        case 311: return "K0";
        case 13: return (pdg > 0) ? "μ+" : "μ-";
        case 11: return (pdg > 0) ? "e+" : "e-";
        case 22: return "γ";
        case 2112: return "Neutron";
        default: return "PDG " + std::to_string(pdg);
    }
}

//********************************************************************
void pdEventDisplay::AddExperimentVariables(OutputManager& output) {
//********************************************************************
    std::cout << "pdEventDisplay::AddExperimentVariables() - Adding ProtoDUNE-specific variables" << std::endl;

    // Use the same tree index as EventDisplayBase
    const Int_t tree_index = 45;

    // Event identification variables
    output.AddVar(tree_index, edrun, "ED_run", "I", "Event Display run number");
    output.AddVar(tree_index, edsubrun, "ED_subrun", "I", "Event Display subrun number");
    output.AddVar(tree_index, edevent, "ED_event", "I", "Event Display event number");
    output.AddVar(tree_index, edhasTrueK0, "ED_hasTrueK0", "I", "Event has true K0");

    // Particle data (variable size arrays)
    output.AddVectorVar(tree_index, edparticle_uniqueID, "ED_particle_uniqueID", "I", "Particle unique IDs", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddVectorVar(tree_index, edparticle_PDG, "ED_particle_PDG", "I", "Particle PDG codes", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddVectorVar(tree_index, edparticle_charge, "ED_particle_charge", "F", "Particle charges", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddVectorVar(tree_index, edparticle_momentum, "ED_particle_momentum", "F", "Particle momenta", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddVectorVar(tree_index, edparticle_nHits, "ED_particle_nHits", "I", "Number of hits per particle", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddVectorVar(tree_index, edparticle_hitIndex, "ED_particle_hitIndex", "I", "Starting hit index for each particle", ednParticles, "ED_nParticles", -kMaxParticles);
    output.AddMatrixVar(tree_index, edparticle_startPos, "ED_particle_startPos", "F", "Particle start positions", ednParticles, "ED_nParticles", -kMaxParticles, 3);
    output.AddMatrixVar(tree_index, edparticle_endPos, "ED_particle_endPos", "F", "Particle end positions", ednParticles, "ED_nParticles", -kMaxParticles, 3);

    // Hit data
    output.AddVectorVar(tree_index, edhit_X, "ED_hit_X", "F", "Hit X positions", edtotalHits, "ED_totalHits", -kMaxHits);
    output.AddVectorVar(tree_index, edhit_Y, "ED_hit_Y", "F", "Hit Y positions", edtotalHits, "ED_totalHits", -kMaxHits);
    output.AddVectorVar(tree_index, edhit_Z, "ED_hit_Z", "F", "Hit Z positions", edtotalHits, "ED_totalHits", -kMaxHits);
    output.AddVectorVar(tree_index, edhit_dEdx, "ED_hit_dEdx", "F", "Hit dEdx values", edtotalHits, "ED_totalHits", -kMaxHits);

    // Call analysis-specific variable addition (e.g., K0 variables in neutralKaonEventDisplay)
    AddAnalysisVariables(output, tree_index);

    std::cout << "pdEventDisplay::AddExperimentVariables() - Variables added successfully" << std::endl;
}

//********************************************************************
void pdEventDisplay::FillExperimentData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
    // Box is cast in analysis-specific derived classes as needed
    (void)box;  // Suppress warning if not used by pdEventDisplay

    // Fill event identification
    output.FillVar(edrun, event.EventInfo->Run);
    output.FillVar(edsubrun, event.EventInfo->SubRun);
    output.FillVar(edevent, event.EventInfo->Event);

    // Check if event has true K0
    Int_t hasTrueK0 = 0;
    for (int i = 0; i < event.nTrueParticles; i++) {
        if (event.TrueParticles[i] && abs(event.TrueParticles[i]->PDG) == 311) {
            hasTrueK0 = 1;
            break;
        }
    }
    output.FillVar(edhasTrueK0, hasTrueK0);

    // Fill particle data
    Int_t hitIndex = 0;

    for (int i = 0; i < event.nParticles && i < kMaxParticles; i++) {
        AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[i]);
        if (!particle) continue;

        output.FillVectorVar(edparticle_uniqueID, particle->UniqueID);

        // Get PDG from TrueObject if available
        Int_t pdg = 0;
        if (particle->TrueObject) {
            AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(particle->TrueObject);
            if (trueParticle) {
                pdg = trueParticle->PDG;
            }
        }
        output.FillVectorVar(edparticle_PDG, (Int_t)pdg);
        output.FillVectorVar(edparticle_charge, (Float_t)particle->Charge);
        output.FillVectorVar(edparticle_momentum, (Float_t)particle->Momentum);

        // Count hits from collection plane (plane 2)
        Int_t nHits = particle->Hits[2].size();
        output.FillVectorVar(edparticle_nHits, (Int_t)nHits);
        output.FillVectorVar(edparticle_hitIndex, (Int_t)hitIndex);

        // Store hit positions for this particle
        for (Int_t h = 0; h < nHits && hitIndex < kMaxHits; h++) {
            const AnaHitPD& hit = particle->Hits[2][nHits - 1 - h]; // Reverse order
            output.FillVectorVar(edhit_X, (Float_t)hit.Position[0]);
            output.FillVectorVar(edhit_Y, (Float_t)hit.Position[1]);
            output.FillVectorVar(edhit_Z, (Float_t)hit.Position[2]);
            output.FillVectorVar(edhit_dEdx, (Float_t)hit.dEdx);
            output.IncrementCounter(edtotalHits);
            hitIndex++;
        }

        // Store positions
        Float_t startPos[3] = {(Float_t)particle->PositionStart[0], (Float_t)particle->PositionStart[1], (Float_t)particle->PositionStart[2]};
        Float_t endPos[3] = {(Float_t)particle->PositionEnd[0], (Float_t)particle->PositionEnd[1], (Float_t)particle->PositionEnd[2]};

        output.FillMatrixVarFromArray(edparticle_startPos, startPos, 3);
        output.FillMatrixVarFromArray(edparticle_endPos, endPos, 3);

        // Increment particle counter
        output.IncrementCounter(ednParticles);
    }

    // Call analysis-specific data filling (e.g., K0 candidates in neutralKaonEventDisplay)
    FillAnalysisData(output, event, box);
}

//********************************************************************
void pdEventDisplay::DrawDetectorGeometry(TEveScene* scene) {
//********************************************************************
    std::cout << "Drawing ProtoDUNE-SP detector geometry..." << std::endl;

    TEveElementList* detectorGeometryGroup = new TEveElementList("Detector Geometry");
    scene->AddElement(detectorGeometryGroup);

    // Create geometry manager if needed (only once per session)
    static bool geoManagerCreated = false;
    if (!gGeoManager && !geoManagerCreated) {
        TGeoManager* geom = new TGeoManager("geom", "ProtoDUNE-SP Detector Geometry");
        geom->SetVisLevel(4);
        geom->SetVisOption(0);
        geoManagerCreated = true;
    }

    // ProtoDUNE-SP dimensions (FINAL CORRECT GEOMETRY)
    // CPA at X=0, APAs at X=-360 and X=+360
    // Positive coordinates only: Y∈[0,600], Z∈[0,696]
    Double_t apaThickness = 2.0;  // Thin in X (drift direction)
    Double_t apaHeight = 600.0;   // Y extent (vertical)
    Double_t apaLength = 696.0;   // Z extent (beam direction)

    Double_t cpaThickness = 2.0;
    Double_t cpaHeight = 600.0;
    Double_t cpaLength = 696.0;

    // Draw CPA (red) - single plane at X=0
    TGeoBBox* cpa = new TGeoBBox("CPA", cpaThickness/2, cpaHeight/2, cpaLength/2);
    TGeoVolume* cpaVol = new TGeoVolume("CPA_Vol", cpa);
    TEveGeoShape* cpaShape = new TEveGeoShape("CPA (X=0)");
    cpaShape->SetShape(cpaVol->GetShape());
    cpaShape->SetMainColor(kRed);
    cpaShape->SetMainTransparency(70);
    cpaShape->RefMainTrans().SetPos(0.0, apaHeight/2, apaLength/2);
    detectorGeometryGroup->AddElement(cpaShape);

    // Draw 6 APAs (blue) - 3 at X=-360 (left side), 3 at X=+360 (right side)
    // Left side APAs (X=-360)
    for (int i = 0; i < 3; i++) {
        Double_t yPos = apaHeight / 2;  // Center at Y=300
        Double_t zPos = i * (apaLength / 3.0) + apaLength / 6.0;  // Distributed along Z

        TGeoBBox* apa = new TGeoBBox(Form("APA_L%d", i), apaThickness/2, apaHeight/2, apaLength/6.0);
        TGeoVolume* apaVol = new TGeoVolume(Form("APA_L%d_Vol", i), apa);
        TEveGeoShape* apaShape = new TEveGeoShape(Form("APA Left %d (X=-360)", i+1));
        apaShape->SetShape(apaVol->GetShape());
        apaShape->SetMainColor(kBlue);
        apaShape->SetMainTransparency(70);
        apaShape->RefMainTrans().SetPos(-360.0, yPos, zPos);
        detectorGeometryGroup->AddElement(apaShape);
    }

    // Right side APAs (X=+360)
    for (int i = 0; i < 3; i++) {
        Double_t yPos = apaHeight / 2;
        Double_t zPos = i * (apaLength / 3.0) + apaLength / 6.0;

        TGeoBBox* apa = new TGeoBBox(Form("APA_R%d", i), apaThickness/2, apaHeight/2, apaLength/6.0);
        TGeoVolume* apaVol = new TGeoVolume(Form("APA_R%d_Vol", i), apa);
        TEveGeoShape* apaShape = new TEveGeoShape(Form("APA Right %d (X=+360)", i+1));
        apaShape->SetShape(apaVol->GetShape());
        apaShape->SetMainColor(kBlue);
        apaShape->SetMainTransparency(70);
        apaShape->RefMainTrans().SetPos(360.0, yPos, zPos);
        detectorGeometryGroup->AddElement(apaShape);
    }

    std::cout << "ProtoDUNE-SP geometry added: 6 APAs, 1 CPA" << std::endl;
}

//********************************************************************
void pdEventDisplay::DrawDetectorGeometry2D(TEveProjectionManager* manager, const std::string& projection_type) {
//********************************************************************
    // 2D projections can show simplified geometry if needed
    // For now, just log that we're skipping 2D geometry
    std::cout << "2D detector geometry for " << projection_type << " (simplified/optional)" << std::endl;
    (void)manager; // Suppress unused warning
}

//********************************************************************
bool pdEventDisplay::ReadEventData(TTree* tree, Int_t run, Int_t subrun, Int_t event) {
//********************************************************************
    // Call base class to find the event
    if (!EventDisplayBase::ReadEventData(tree, run, subrun, event)) {
        return false;
    }

    // Set branch addresses for ProtoDUNE-specific data
    tree->SetBranchAddress("ED_hasTrueK0", &_hasTrueK0);
    tree->SetBranchAddress("ED_nParticles", &_nParticles);
    tree->SetBranchAddress("ED_totalHits", &_totalHits);
    tree->SetBranchAddress("ED_nK0Candidates", &_nK0Candidates);

    tree->SetBranchAddress("ED_particle_uniqueID", _particle_uniqueID);
    tree->SetBranchAddress("ED_particle_PDG", _particle_PDG);
    tree->SetBranchAddress("ED_particle_charge", _particle_charge);
    tree->SetBranchAddress("ED_particle_momentum", _particle_momentum);
    tree->SetBranchAddress("ED_particle_nHits", _particle_nHits);
    tree->SetBranchAddress("ED_particle_hitIndex", _particle_hitIndex);
    tree->SetBranchAddress("ED_particle_startPos", _particle_startPos);
    tree->SetBranchAddress("ED_particle_endPos", _particle_endPos);

    tree->SetBranchAddress("ED_hit_X", _hit_X);
    tree->SetBranchAddress("ED_hit_Y", _hit_Y);
    tree->SetBranchAddress("ED_hit_Z", _hit_Z);
    tree->SetBranchAddress("ED_hit_dEdx", _hit_dEdx);

    tree->SetBranchAddress("ED_k0_creationVtxPos", _k0_creationVtxPos);
    tree->SetBranchAddress("ED_k0_annihilationVtxPos", _k0_annihilationVtxPos);
    tree->SetBranchAddress("ED_k0_startPos", _k0_startPos);
    tree->SetBranchAddress("ED_k0_endPos", _k0_endPos);
    tree->SetBranchAddress("ED_k0_daughter1ID", _k0_daughter1ID);
    tree->SetBranchAddress("ED_k0_daughter2ID", _k0_daughter2ID);
    tree->SetBranchAddress("ED_k0_parentID", _k0_parentID);
    tree->SetBranchAddress("ED_k0_creationVtxRadius", _k0_creationVtxRadius);
    tree->SetBranchAddress("ED_k0_annihilationVtxRadius", _k0_annihilationVtxRadius);

    // Read the current entry
    tree->GetEntry(tree->GetReadEntry());

    std::cout << "Read event data: " << _nParticles << " particles, "
              << _totalHits << " hits, " << _nK0Candidates << " K0 candidates" << std::endl;

    // Read analysis-specific data
    if (!ReadAnalysisData(tree)) {
        std::cerr << "Failed to read analysis-specific data" << std::endl;
        return false;
    }

    return true;
}

//********************************************************************
void pdEventDisplay::DrawParticles3D(TEveScene* scene) {
//********************************************************************
    std::cout << "Drawing particles in 3D..." << std::endl;

    TEveElementList* particlesGroup = new TEveElementList("Particles");
    scene->AddElement(particlesGroup);

    auto isDuplicateK0RecoTrajectory = [&](Int_t particleIndex, Int_t pdg, Int_t nHits) {
        if (nHits > 0) return false;
        const Int_t absPdg = std::abs(pdg);
        if (!(absPdg == 311 || absPdg == 310 || absPdg == 130)) return false;

        const Float_t startX = _particle_startPos[particleIndex][0];
        const Float_t startY = _particle_startPos[particleIndex][1];
        const Float_t startZ = _particle_startPos[particleIndex][2];
        const Float_t endX = _particle_endPos[particleIndex][0];
        const Float_t endY = _particle_endPos[particleIndex][1];
        const Float_t endZ = _particle_endPos[particleIndex][2];
        if (startX <= -900 || startY <= -900 || startZ <= -900 ||
            endX <= -900 || endY <= -900 || endZ <= -900) {
            return false;
        }

        const Float_t tol2 = 1e-4f;
        auto samePoint = [&](Float_t ax, Float_t ay, Float_t az, Float_t bx, Float_t by, Float_t bz) {
            const Float_t dx = ax - bx;
            const Float_t dy = ay - by;
            const Float_t dz = az - bz;
            return (dx * dx + dy * dy + dz * dz) < tol2;
        };

        for (Int_t k = 0; k < _nK0Candidates && k < kMaxK0; ++k) {
            if (_k0_startPos[k][0] <= -900 || _k0_endPos[k][0] <= -900) continue;
            if (samePoint(startX, startY, startZ,
                          _k0_startPos[k][0], _k0_startPos[k][1], _k0_startPos[k][2]) &&
                samePoint(endX, endY, endZ,
                          _k0_endPos[k][0], _k0_endPos[k][1], _k0_endPos[k][2])) {
                return true;
            }
        }
        return false;
    };

    Int_t hitIdx = 0;
    for (Int_t i = 0; i < _nParticles && i < kMaxParticles; i++) {
        Int_t nHits = _particle_nHits[i];
        Int_t pdg = _particle_PDG[i];
        Int_t uniqueID = _particle_uniqueID[i];

        if (isDuplicateK0RecoTrajectory(i, pdg, nHits)) {
            hitIdx += nHits;
            continue;
        }

        Int_t color = GetParticleColor(pdg);
        std::string pdgName = GetParticleName(pdg);

        TEveElementList* particleGroup =
            new TEveElementList(Form("Particle UID=%d TruePDG=%d (%s)", uniqueID, pdg, pdgName.c_str()));
        particlesGroup->AddElement(particleGroup);

        std::map<Int_t, TEvePointSet*> hitsBySizeBin;
        TEvePointSet* firstHitSet = nullptr;
        TEvePointSet* startPosSet = nullptr;
        TEvePointSet* endPosSet = nullptr;

        // Check if particle has hits
        if (nHits > 0) {
            // Particle has hits - draw hit-by-hit trajectory

            firstHitSet = new TEvePointSet("First Hit");
            firstHitSet->SetMarkerStyle(29);
            firstHitSet->SetMarkerSize(2.0);
            firstHitSet->SetMainColor(color);
            particleGroup->AddElement(firstHitSet);

            // Add this particle's hits
            for (Int_t h = 0; h < nHits; h++) {
                if (hitIdx < _totalHits) {
                    Float_t x = _hit_X[hitIdx];
                    Float_t y = _hit_Y[hitIdx];
                    Float_t z = _hit_Z[hitIdx];
                    Float_t dEdx = _hit_dEdx[hitIdx];

                    if (x > -900 && y > -900 && z > -900) {
                        if (h == nHits - 1) {
                            // Last hit in array = first hit spatially
                            firstHitSet->SetNextPoint(x, y, z);
                        } else {
                            const Int_t sizeBin = HitSizeBinFromdEdx(dEdx);
                            if (hitsBySizeBin.find(sizeBin) == hitsBySizeBin.end()) {
                                TEvePointSet* hitSet =
                                    new TEvePointSet(Form("Hits dEdxBin%d", sizeBin));
                                hitSet->SetMarkerStyle(20);
                                hitSet->SetMarkerSize(HitSizeForBin(sizeBin));
                                hitSet->SetMainColor(color);
                                hitsBySizeBin[sizeBin] = hitSet;
                                particleGroup->AddElement(hitSet);
                            }
                            hitsBySizeBin[sizeBin]->SetNextPoint(x, y, z);
                        }
                    }
                    hitIdx++;
                }
            }

            startPosSet = new TEvePointSet("Pandora Start");
            startPosSet->SetMarkerStyle(21);
            startPosSet->SetMarkerSize(1.2);
            startPosSet->SetMainColor(color);
            particleGroup->AddElement(startPosSet);

            endPosSet = new TEvePointSet("Pandora End");
            endPosSet->SetMarkerStyle(22);
            endPosSet->SetMarkerSize(1.2);
            endPosSet->SetMainColor(color);
            particleGroup->AddElement(endPosSet);
        } else {
            // Particle has no hits - draw simple line from start to end position
            if (_particle_startPos[i][0] > -900 && _particle_endPos[i][0] > -900) {
                TEveLine* line = new TEveLine("Line");
                line->SetNextPoint(_particle_startPos[i][0],
                                  _particle_startPos[i][1],
                                  _particle_startPos[i][2]);
                line->SetNextPoint(_particle_endPos[i][0],
                                  _particle_endPos[i][1],
                                  _particle_endPos[i][2]);
                line->SetLineColor(color);
                line->SetLineWidth(2);
                line->SetLineStyle(2); // Dashed for no-hit particles
                particleGroup->AddElement(line);
            }

            startPosSet = new TEvePointSet("Pandora Start");
            startPosSet->SetMarkerStyle(21);
            startPosSet->SetMarkerSize(1.2);
            startPosSet->SetMainColor(color);
            particleGroup->AddElement(startPosSet);

            endPosSet = new TEvePointSet("Pandora End");
            endPosSet->SetMarkerStyle(22);
            endPosSet->SetMarkerSize(1.2);
            endPosSet->SetMainColor(color);
            particleGroup->AddElement(endPosSet);
        }

        if (startPosSet &&
            _particle_startPos[i][0] > -900 && _particle_startPos[i][1] > -900 && _particle_startPos[i][2] > -900) {
            startPosSet->SetNextPoint(_particle_startPos[i][0],
                                      _particle_startPos[i][1],
                                      _particle_startPos[i][2]);
        }

        if (endPosSet &&
            _particle_endPos[i][0] > -900 && _particle_endPos[i][1] > -900 && _particle_endPos[i][2] > -900) {
            endPosSet->SetNextPoint(_particle_endPos[i][0],
                                    _particle_endPos[i][1],
                                    _particle_endPos[i][2]);
        }
    }

    std::cout << "3D particles drawn" << std::endl;

    // Draw analysis-specific content
    DrawAnalysisContent3D(scene);
}

//********************************************************************
void pdEventDisplay::DrawParticles2D(TEveProjectionManager* manager, const std::string& projection_type) {
//********************************************************************
    std::cout << "Drawing particles in 2D (" << projection_type << ")..." << std::endl;

    // Project elements from 3D scene
    // Note: This is a simplified implementation
    // A full implementation would create 2D-specific elements

    (void)manager; // Suppress unused warning for now
    (void)projection_type;

    std::cout << "2D particle projection (using TEve auto-projection)" << std::endl;

    // Draw analysis-specific content
    DrawAnalysisContent2D(manager, projection_type);
}

// ========== 2D TCanvas Drawing Methods ==========

//********************************************************************
void pdEventDisplay::DrawDetectorCanvas2D(TCanvas* canvas, const std::string& projection_type) {
//********************************************************************
    canvas->cd();

    Float_t xmin, xmax, ymin, ymax;
    std::string xtitle, ytitle;

    if (projection_type == "XY") {
        xmin = -400; xmax = 400;
        ymin = -50;  ymax = 650;
        xtitle = "X (cm - drift direction)";
        ytitle = "Y (cm - vertical)";
    } else if (projection_type == "XZ") {
        xmin = -400; xmax = 400;
        ymin = -50;  ymax = 750;
        xtitle = "X (cm - drift direction)";
        ytitle = "Z (cm - beam direction)";
    } else if (projection_type == "YZ") {
        xmin = -50;  xmax = 650;
        ymin = -50;  ymax = 750;
        xtitle = "Y (cm - vertical)";
        ytitle = "Z (cm - beam direction)";
    } else {
        return;
    }

    TH2F* h2 = new TH2F(Form("h2_base_%s", projection_type.c_str()),
                         Form("ProtoDUNE - %s View;%s;%s",
                              projection_type.c_str(), xtitle.c_str(), ytitle.c_str()),
                         100, xmin, xmax, 100, ymin, ymax);
    h2->SetStats(false);
    h2->SetDirectory(0);
    h2->GetXaxis()->SetLabelSize(0.04);
    h2->GetYaxis()->SetLabelSize(0.04);
    h2->GetXaxis()->SetTitleSize(0.045);
    h2->GetYaxis()->SetTitleSize(0.045);
    h2->Draw();

    const Float_t apaH = 600.0, apaL = 696.0;

    if (projection_type == "XY") {
        TLine* cpa = new TLine(0, 0, 0, apaH);
        cpa->SetLineColor(kRed); cpa->SetLineWidth(3); cpa->Draw("SAME");
        TLine* apa1 = new TLine(-360, 0, -360, apaH);
        apa1->SetLineColor(kBlue); apa1->SetLineWidth(3); apa1->Draw("SAME");
        TLine* apa2 = new TLine(360, 0, 360, apaH);
        apa2->SetLineColor(kBlue); apa2->SetLineWidth(3); apa2->Draw("SAME");
    } else if (projection_type == "XZ") {
        TLine* cpa = new TLine(0, 0, 0, apaL);
        cpa->SetLineColor(kRed); cpa->SetLineWidth(3); cpa->Draw("SAME");
        TLine* apa1 = new TLine(-360, 0, -360, apaL);
        apa1->SetLineColor(kBlue); apa1->SetLineWidth(3); apa1->Draw("SAME");
        TLine* apa2 = new TLine(360, 0, 360, apaL);
        apa2->SetLineColor(kBlue); apa2->SetLineWidth(3); apa2->Draw("SAME");
    } else if (projection_type == "YZ") {
        TLine* bottom = new TLine(0, 0, 0, apaL);
        TLine* top    = new TLine(apaH, 0, apaH, apaL);
        TLine* left   = new TLine(0, 0, apaH, 0);
        TLine* right  = new TLine(0, apaL, apaH, apaL);
        for (auto* l : {bottom, top, left, right}) {
            l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw("SAME");
        }
    }
}

//********************************************************************
void pdEventDisplay::DrawParticlesCanvas2D(TCanvas* canvas, const std::string& projection_type) {
//********************************************************************
    canvas->cd();

    auto project = [&](Float_t x, Float_t y, Float_t z,
                       Float_t& c1, Float_t& c2) -> bool {
        if (projection_type == "XY") {
            if (x <= -900 || y <= -900) return false;
            c1 = x; c2 = y; return true;
        } else if (projection_type == "XZ") {
            if (x <= -900 || z <= -900) return false;
            c1 = x; c2 = z; return true;
        } else if (projection_type == "YZ") {
            if (y <= -900 || z <= -900) return false;
            c1 = y; c2 = z; return true;
        }
        return false;
    };

    std::set<Int_t> uniquePDGs;
    std::map<Int_t, TGraph*> graphsByPDG;

    Int_t hitIdx = 0;
    for (Int_t i = 0; i < _nParticles && i < kMaxParticles; i++) {
        Int_t nHits = _particle_nHits[i];
        Int_t pdg   = _particle_PDG[i];
        Int_t color = GetParticleColor(pdg);
        uniquePDGs.insert(pdg);

        if (nHits > 0) {
            for (Int_t h = 0; h < nHits; h++) {
                if (hitIdx >= _totalHits) break;
                Float_t c1, c2;
                if (!project(_hit_X[hitIdx], _hit_Y[hitIdx], _hit_Z[hitIdx], c1, c2)) {
                    hitIdx++; continue;
                }
                const Float_t dEdx = _hit_dEdx[hitIdx];

                Int_t style = (h == 0) ? 29 : 20;
                Float_t size = (h == 0) ? 1.2f : HitMarkerSizeFromdEdx(dEdx);
                TMarker* m = new TMarker(c1, c2, style);
                m->SetMarkerColor(color);
                m->SetMarkerSize(size);
                m->Draw("SAME");
                hitIdx++;
            }

            if (graphsByPDG.find(pdg) == graphsByPDG.end()) {
                TGraph* gr = new TGraph(1);
                gr->SetPoint(0, -9999, -9999);
                gr->SetMarkerStyle(20); gr->SetMarkerSize(0.8);
                gr->SetMarkerColor(color);
                graphsByPDG[pdg] = gr;
            }
        } else {
            Float_t x1, y1, x2, y2;
            bool s = project(_particle_startPos[i][0], _particle_startPos[i][1],
                             _particle_startPos[i][2], x1, y1);
            bool e = project(_particle_endPos[i][0], _particle_endPos[i][1],
                             _particle_endPos[i][2], x2, y2);
            if (s && e) {
                TLine* line = new TLine(x1, y1, x2, y2);
                line->SetLineColor(color); line->SetLineStyle(2);
                line->SetLineWidth(2); line->Draw("SAME");
                TMarker* sm = new TMarker(x1, y1, 29);
                sm->SetMarkerColor(color); sm->SetMarkerSize(1.2f);
                sm->Draw("SAME");
            }
        }
    }

    TLegend* legend = new TLegend(0.85, 0.7, 0.98, 0.95);
    legend->SetTextSize(0.025);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(1);
    for (Int_t pdg : uniquePDGs) {
        std::string name = GetParticleName(pdg);
        Int_t color = GetParticleColor(pdg);
        if (graphsByPDG.count(pdg)) {
            legend->AddEntry(graphsByPDG[pdg], name.c_str(), "p");
        } else {
            TMarker* dm = new TMarker(0, 0, 20);
            dm->SetMarkerColor(color);
            legend->AddEntry(dm, name.c_str(), "p");
        }
    }
    legend->Draw("SAME");
}

// ========== Analysis-Specific Virtual Methods (Default Implementations) ==========

//********************************************************************
void pdEventDisplay::AddAnalysisVariables(OutputManager& output, Int_t tree_index) {
//********************************************************************
    // Default: no analysis-specific variables
    // Override in derived classes (e.g., neutralKaonEventDisplay)
    (void)output;
    (void)tree_index;
}

//********************************************************************
void pdEventDisplay::FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
    // Default: no analysis-specific data
    (void)output;
    (void)event;
    (void)box;
}

//********************************************************************
bool pdEventDisplay::ReadAnalysisData(TTree* tree) {
//********************************************************************
    // Default: no analysis-specific data to read
    (void)tree;
    return true;
}

//********************************************************************
void pdEventDisplay::DrawAnalysisContent3D(TEveScene* scene) {
//********************************************************************
    // Default: no analysis-specific content
    (void)scene;
}

//********************************************************************
void pdEventDisplay::DrawAnalysisContent2D(TEveProjectionManager* manager, const std::string& projection_type) {
//********************************************************************
    // Default: no analysis-specific content
    (void)manager;
    (void)projection_type;
}
