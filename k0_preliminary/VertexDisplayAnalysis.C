// File: VertexDisplayAnalysis.C
// Combines event iteration structure from PreliminaryK0Selection with vertex creation algorithm from neutralKaonAnalysis
// Displays events and draws the vertices found in each event

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TBox.h"
#include "TMarker.h"
#include "TText.h"
#include "TLegend.h"
#include "TStyle.h"

// HIGHLAND includes - these will be available when libraries are loaded
// #include "AnaSpill.h"
// #include "AnaBeam.h"
// #include "AnaBunch.h"
// #include "AnaParticlePD.h"
// #include "AnaTrueParticleB.h"
// #include "AnaVertexPD.h"
// #include "AnaTrueVertexPD.h"

// Vertex structure for display
struct DisplayVertex {
    double position[3];
    std::vector<AnaParticlePD*> particles;
    AnaParticlePD* parent;
    int nParticles;
    bool isTrueVertex;
    int vertexType; // 0: reconstructed, 1: true
};

// Forward declarations
std::vector<DisplayVertex> CreateVertices(std::vector<AnaParticlePD*>& validParticles, AnaTrueParticleB* beamTrue);
void CreateEventVisualization(const std::vector<DisplayVertex>& vertices, const std::vector<AnaParticlePD*>& allParticles, AnaTrueParticleB* beamTrue, int eventNumber, const std::vector<Int_t>& k0ParentIDs, AnaSpill* spill);

void VertexDisplayAnalysis(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA/6GeV_prod4a_01_minitree_2023-01-27.root", int maxEvents = 1000) {

    std::cout << "=== Vertex Display Analysis ===" << std::endl;
    std::cout << "Input file: " << filename << std::endl;
    std::cout << "Maximum events to process: " << maxEvents << std::endl;

    // Open the file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // Get the tree
    TTree *MiniTree = (TTree*)f->Get("MiniTree");
    if (!MiniTree) {
        std::cerr << "Could not find MiniTree in file" << std::endl;
        return;
    }

    // Set up the object for the branch
    AnaSpill* spill = nullptr;
    MiniTree->SetBranchAddress("Spill", &spill);

    Long64_t nEntries = MiniTree->GetEntries();
    std::cout << "Total events in file: " << nEntries << std::endl;

    int eventsProcessed = 0;

    for (Long64_t i = 0; i < nEntries && eventsProcessed < maxEvents; i++) {
        MiniTree->GetEntry(i);

        if (!spill) continue;

        std::cout << "\n=== Processing Event " << i << " ===" << std::endl;

        // Get beam particle
        AnaParticleB* beamParticle = nullptr;
        AnaTrueParticleB* beamTrue = nullptr;
        if (spill->Beam) {
            AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
            if (beam && beam->BeamParticle) {
                beamParticle = beam->BeamParticle;
                beamTrue = beamParticle->GetTrueParticle();
            }
        }

        if (!beamTrue) {
            std::cout << "Skipping event - no beam particle or true particle" << std::endl;
            continue;
        }

        std::cout << "Beam particle: PDG=" << beamTrue->PDG << ", ID=" << beamTrue->ID << std::endl;
        std::cout << "Beam end position: (" << beamTrue->PositionEnd[0] << ", " << beamTrue->PositionEnd[1] << ", " << beamTrue->PositionEnd[2] << ")" << std::endl;

        // Collect all valid reconstructed particles
        std::vector<AnaParticlePD*> validParticles;

        for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
            AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
            if (!bunch) continue;

            for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                AnaParticleB* particle = bunch->Particles[j];
                if (!particle) continue;

                AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                if (!particlePD) continue;

                // Check if particle has valid start and end positions
                if (particlePD->PositionStart[0] < -900 || particlePD->PositionStart[1] < -900 || particlePD->PositionStart[2] < -900 ||
                    particlePD->PositionEnd[0] < -900 || particlePD->PositionEnd[1] < -900 || particlePD->PositionEnd[2] < -900) {
                    continue;
                }

                validParticles.push_back(particlePD);
            }
        }

        std::cout << "Found " << validParticles.size() << " valid reconstructed particles" << std::endl;

        if (validParticles.size() < 2) {
            std::cout << "Not enough particles to form vertices" << std::endl;
            continue;
        }

        // Apply vertex creation algorithm from neutralKaonAnalysis
        std::vector<DisplayVertex> foundVertices = CreateVertices(validParticles, beamTrue);

        std::cout << "Found " << foundVertices.size() << " vertices in this event" << std::endl;

        // First loop: Group particles into vertices and identify K0 ParentIDs
        bool shouldDisplayEvent = false;
        std::vector<Int_t> k0ParentIDs; // Store ParentIDs of K0 particles

        for (const auto& vertex : foundVertices) {
            // Check if this vertex has a K+ parent (PDG=321)
            AnaTrueParticleB* parentTrue = vertex.parent->GetTrueParticle();
            if (!parentTrue || parentTrue->PDG != 321) continue; // Skip if not K+

            // Count pions in this vertex
            int nPions = 0;
            std::vector<AnaTrueParticleB*> pionParticles;

            for (size_t j = 1; j < vertex.particles.size(); j++) { // Skip parent (index 0)
                AnaTrueParticleB* daughterTrue = vertex.particles[j]->GetTrueParticle();
                if (!daughterTrue) continue;

                // Check if this is a pion (π+ = 211, π- = -211)
                if (daughterTrue->PDG == 211 || daughterTrue->PDG == -211) {
                    nPions++;
                    pionParticles.push_back(daughterTrue);
                }
            }

            // Check if we have at least 2 pions
            if (nPions >= 2) {
                // Check if at least 2 pions come from a K0 (PDG=310)
                int pionsFromK0 = 0;
                std::vector<AnaTrueParticleB*> k0Pions;

                for (const auto& pion : pionParticles) {
                    // Check if this pion's parent is a K0 (PDG=310) using true particle information
                    if (pion->ParentPDG == 310) {
                        pionsFromK0++;
                        k0Pions.push_back(pion);

                        // Store the K0 ParentID for later lookup
                        bool alreadyStored = false;
                        for (const auto& existingID : k0ParentIDs) {
                            if (existingID == pion->ParentID) {
                                alreadyStored = true;
                                break;
                            }
                        }
                        if (!alreadyStored) {
                            k0ParentIDs.push_back(pion->ParentID);
                        }
                    }
                }

                if (pionsFromK0 >= 2) {
                    shouldDisplayEvent = true;
                    std::cout << "*** INTERESTING EVENT FOUND ***" << std::endl;
                    std::cout << "  K+ vertex with K0→π+π- decay detected!" << std::endl;
                    std::cout << "  Parent K+: PDG=" << parentTrue->PDG << ", ID=" << parentTrue->ID << std::endl;
                    std::cout << "  Found " << pionsFromK0 << " pions from K0 out of " << nPions << " total pions" << std::endl;
                    for (size_t k = 0; k < k0Pions.size(); k++) {
                        std::cout << "  K0 Daughter π" << (k0Pions[k]->PDG > 0 ? "+" : "-")
                                  << ": PDG=" << k0Pions[k]->PDG
                                  << ", ID=" << k0Pions[k]->ID
                                  << ", ParentID=" << k0Pions[k]->ParentID << std::endl;
                    }
                    break; // Found what we're looking for, no need to check other vertices
                }
            }
        }

        // Skip this event if it doesn't meet our criteria
        if (!shouldDisplayEvent) {
            std::cout << "Event " << i << " does not contain K+→K0→π+π- decay, skipping visualization." << std::endl;
            continue;
        }

        std::cout << "*** DISPLAYING INTERESTING EVENT ***" << std::endl;

        // Display vertex information
        for (size_t v = 0; v < foundVertices.size(); v++) {
            const DisplayVertex& vertex = foundVertices[v];
            std::cout << "  Vertex " << v << ": position=(" << vertex.position[0] << ", " << vertex.position[1] << ", " << vertex.position[2] << "), nParticles=" << vertex.nParticles << std::endl;

            // Show particles in this vertex (TRUE PDG for identification, RECO positions for vertex creation)
            for (size_t p = 0; p < vertex.particles.size(); p++) {
                AnaParticlePD* part = vertex.particles[p];
                AnaTrueParticleB* truePart = part->GetTrueParticle();
                if (truePart) {
                    std::cout << "    Particle " << p << ": PDG=" << truePart->PDG << ", TrueID=" << truePart->ID;
                    if (p == 0) std::cout << " (parent)";
                    std::cout << std::endl;
                }
            }
        }

        // Create visualization
        if (foundVertices.size() > 0) {
            try {
                CreateEventVisualization(foundVertices, validParticles, beamTrue, i, k0ParentIDs, spill);
            } catch (const std::exception& e) {
                std::cerr << "Error creating visualization for event " << i << ": " << e.what() << std::endl;
                continue;
            }
        }

        eventsProcessed++;
    }

    std::cout << "\nProcessed " << eventsProcessed << " events with vertices" << std::endl;
    f->Close();
    delete f;
}

// Vertex creation algorithm based on neutralKaonAnalysis
std::vector<DisplayVertex> CreateVertices(std::vector<AnaParticlePD*>& validParticles, AnaTrueParticleB* beamTrue) {
    std::vector<DisplayVertex> vertices;

    const double maxVertexRadius = 30.0; // cm - maximum radius for a vertex
    const int minParticlesPerVertex = 2; // minimum particles required to form a vertex

    // Track which particles have been assigned to vertices
    std::vector<bool> particleAssigned(validParticles.size(), false);

    // First, find all particles with valid end positions (potential parents)
    std::vector<AnaParticlePD*> potentialParents;
    for (size_t i = 0; i < validParticles.size(); i++) {
        AnaParticlePD* part = validParticles[i];
        if (!part) continue;

        // Check if this particle has valid end position
        if (part->PositionEnd[0] < -900 || part->PositionEnd[1] < -900 || part->PositionEnd[2] < -900) {
            continue;
        }

        potentialParents.push_back(part);
    }

    std::cout << "Found " << potentialParents.size() << " potential parent particles" << std::endl;

    // For each potential parent, find daughters whose start positions are within maxVertexRadius of parent's end position
    for (const auto& parent : potentialParents) {
        std::vector<AnaParticlePD*> vertexParticles;

        // Add the parent to the vertex (it will be the parent particle)
        vertexParticles.push_back(parent);

        // Find all particles whose start positions are within maxVertexRadius of parent's end position
        for (size_t p = 0; p < validParticles.size(); p++) {
            if (particleAssigned[p]) continue; // Skip already assigned particles

            AnaParticlePD* candidateParticle = validParticles[p];
            if (!candidateParticle) continue;

            // Skip if this is the same particle as the parent
            if (candidateParticle == parent) continue;

            // Calculate distance from parent's end position to candidate's start position
            double dist = sqrt((candidateParticle->PositionStart[0] - parent->PositionEnd[0]) * (candidateParticle->PositionStart[0] - parent->PositionEnd[0]) +
                              (candidateParticle->PositionStart[1] - parent->PositionEnd[1]) * (candidateParticle->PositionStart[1] - parent->PositionEnd[1]) +
                              (candidateParticle->PositionStart[2] - parent->PositionEnd[2]) * (candidateParticle->PositionStart[2] - parent->PositionEnd[2]));

            if (dist <= maxVertexRadius) {
                vertexParticles.push_back(candidateParticle);
                particleAssigned[p] = true;
            }
        }

        // Only create a vertex if we have enough particles (parent + daughters)
        if (vertexParticles.size() >= minParticlesPerVertex) {
            DisplayVertex vertex;

            // Set vertex position to the parent's end position
            vertex.position[0] = parent->PositionEnd[0];
            vertex.position[1] = parent->PositionEnd[1];
            vertex.position[2] = parent->PositionEnd[2];

            vertex.particles = vertexParticles;
            vertex.parent = parent;
            vertex.nParticles = vertexParticles.size();
            vertex.isTrueVertex = false; // For now, we'll mark as reconstructed
            vertex.vertexType = 0; // Reconstructed vertex

            vertices.push_back(vertex);
        }
    }

    return vertices;
}

// Create visualization of the event (following PreliminaryK0Selection pattern)
void CreateEventVisualization(const std::vector<DisplayVertex>& vertices, const std::vector<AnaParticlePD*>& allParticles, AnaTrueParticleB* beamTrue, int eventNumber, const std::vector<Int_t>& k0ParentIDs, AnaSpill* spill) {

    // Create canvas title with vertex information
    std::string canvasTitle = "Event " + std::to_string(eventNumber) + " - Vertex Display Analysis";
    canvasTitle += " - Vertices: " + std::to_string(vertices.size()) + " (30cm radius), Particles: " + std::to_string(allParticles.size());

    // Create unique canvas name for each event
    std::string canvasName = "c1_event_" + std::to_string(eventNumber);
    TCanvas *c1 = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1400, 700);
    c1->Divide(2, 1); // Two subplots side by side (XY and XZ projections)

    // Define particle colors (similar to PreliminaryK0Selection)
    std::map<int, int> particleColors;
    particleColors[211] = kBlue;    // π+
    particleColors[-211] = kOrange; // π-
    particleColors[111] = kCyan;    // π0
    particleColors[321] = kRed;     // K+
    particleColors[-321] = kTeal;   // K-
    particleColors[311] = kYellow+2; // K0
    particleColors[-311] = kTeal+2;  // K̄0
    particleColors[310] = kMagenta;  // K0S
    particleColors[13] = kBlack;     // μ-
    particleColors[-13] = kViolet+2; // μ+
    particleColors[11] = kPink;      // e-
    particleColors[-11] = kPink+2;   // e+
    particleColors[2212] = kGray;    // p
    particleColors[-2212] = kGray+1; // p̄
    particleColors[2112] = kGray+2;  // n
    particleColors[-2112] = kGray+3; // n̄

    // Create new graphs for this event only (following PreliminaryK0Selection pattern)
    std::map<int, TGraph*> eventParticleGraphs;

    // Process reconstructed particles to create initial graphs (using TRUE PDG for identification, RECO positions for display)
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue; // Skip particles without true information

        // Use TRUE particle PDG for identification and coloring
        int particleType = truePart->PDG;

        // Create graph for this particle type if it doesn't exist
        if (eventParticleGraphs.find(particleType) == eventParticleGraphs.end()) {
            eventParticleGraphs[particleType] = new TGraph();
            eventParticleGraphs[particleType]->SetMarkerStyle(20);
            eventParticleGraphs[particleType]->SetMarkerSize(0.5);
            if (particleColors.find(particleType) != particleColors.end()) {
                eventParticleGraphs[particleType]->SetMarkerColor(particleColors[particleType]);
            } else {
                eventParticleGraphs[particleType]->SetMarkerColor(kGray);
            }
        }

        // Add hits from all planes (using RECONSTRUCTED hit positions)
        for (int plane = 0; plane < 3; plane++) {
            for (UInt_t h = 0; h < particle->Hits[plane].size(); h++) {
                AnaHitPD hit = particle->Hits[plane][h];
                if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                    int nPoints = eventParticleGraphs[particleType]->GetN();
                    eventParticleGraphs[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Y());
                }
            }
        }
    }

    // Create graphs for both projections (following PreliminaryK0Selection pattern)
    std::map<int, TGraph*> eventParticleGraphsXY;
    std::map<int, TGraph*> eventParticleGraphsXZ;

    // Copy graphs for XY projection
    for (auto& pair : eventParticleGraphs) {
        eventParticleGraphsXY[pair.first] = new TGraph();
        eventParticleGraphsXY[pair.first]->SetMarkerStyle(20);
        eventParticleGraphsXY[pair.first]->SetMarkerSize(0.5);
        eventParticleGraphsXY[pair.first]->SetMarkerColor(pair.second->GetMarkerColor());
    }

    // Copy graphs for XZ projection
    for (auto& pair : eventParticleGraphs) {
        eventParticleGraphsXZ[pair.first] = new TGraph();
        eventParticleGraphsXZ[pair.first]->SetMarkerStyle(20);
        eventParticleGraphsXZ[pair.first]->SetMarkerSize(0.5);
        eventParticleGraphsXZ[pair.first]->SetMarkerColor(pair.second->GetMarkerColor());
    }

    // Special graphs for vertices
    TGraph* vertexGraphXY = new TGraph();
    vertexGraphXY->SetMarkerStyle(21);
    vertexGraphXY->SetMarkerSize(1.5);
    vertexGraphXY->SetMarkerColor(kBlue);

    TGraph* vertexGraphXZ = new TGraph();
    vertexGraphXZ->SetMarkerStyle(21);
    vertexGraphXZ->SetMarkerSize(1.5);
    vertexGraphXZ->SetMarkerColor(kBlue);

    // Special graphs for parent particle end positions (reconstructed end position)
    TGraph* parentEndGraphXY = new TGraph();
    parentEndGraphXY->SetMarkerStyle(22);
    parentEndGraphXY->SetMarkerSize(1.8);
    parentEndGraphXY->SetMarkerColor(kRed);

    TGraph* parentEndGraphXZ = new TGraph();
    parentEndGraphXZ->SetMarkerStyle(22);
    parentEndGraphXZ->SetMarkerSize(1.8);
    parentEndGraphXZ->SetMarkerColor(kRed);

    // Special graphs for parent particle last hits (actual last hit position)
    TGraph* parentLastHitGraphXY = new TGraph();
    parentLastHitGraphXY->SetMarkerStyle(23);
    parentLastHitGraphXY->SetMarkerSize(1.8);
    parentLastHitGraphXY->SetMarkerColor(kOrange);

    TGraph* parentLastHitGraphXZ = new TGraph();
    parentLastHitGraphXZ->SetMarkerStyle(23);
    parentLastHitGraphXZ->SetMarkerSize(1.8);
    parentLastHitGraphXZ->SetMarkerColor(kOrange);

    // Special graphs for daughter particle first hits (actual first hit position)
    TGraph* daughterFirstHitGraphXY = new TGraph();
    daughterFirstHitGraphXY->SetMarkerStyle(24);
    daughterFirstHitGraphXY->SetMarkerSize(1.5);
    daughterFirstHitGraphXY->SetMarkerColor(kGreen);

    TGraph* daughterFirstHitGraphXZ = new TGraph();
    daughterFirstHitGraphXZ->SetMarkerStyle(24);
    daughterFirstHitGraphXZ->SetMarkerSize(1.5);
    daughterFirstHitGraphXZ->SetMarkerColor(kGreen);

    // Special graphs for daughter particle start positions (reconstructed start position)
    TGraph* daughterStartGraphXY = new TGraph();
    daughterStartGraphXY->SetMarkerStyle(25);
    daughterStartGraphXY->SetMarkerSize(1.5);
    daughterStartGraphXY->SetMarkerColor(kMagenta);

    TGraph* daughterStartGraphXZ = new TGraph();
    daughterStartGraphXZ->SetMarkerStyle(25);
    daughterStartGraphXZ->SetMarkerSize(1.5);
    daughterStartGraphXZ->SetMarkerColor(kMagenta);

    // Store vertex circles for cleanup
    std::vector<TEllipse*> vertexCirclesXY;
    std::vector<TEllipse*> vertexCirclesXZ;

    // Store particle labels for cleanup
    std::vector<TText*> particleLabelsXY;
    std::vector<TText*> particleLabelsXZ;

    // Store K0 trajectories for cleanup
    std::vector<TLine*> k0TrajectoriesXY;
    std::vector<TLine*> k0TrajectoriesXZ;


    // Fill XY and XZ graphs with hit data (using TRUE PDG for identification, RECO positions for display)
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue; // Skip particles without true information

        // Use TRUE particle PDG for identification
        int particleType = truePart->PDG;

        // Add hits to XY projection (using RECONSTRUCTED hit positions)
        if (eventParticleGraphsXY.find(particleType) != eventParticleGraphsXY.end()) {
            for (int plane = 0; plane < 3; plane++) {
                for (UInt_t h = 0; h < particle->Hits[plane].size(); h++) {
                    AnaHitPD hit = particle->Hits[plane][h];
                    if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                        int nPointsXY = eventParticleGraphsXY[particleType]->GetN();
                        eventParticleGraphsXY[particleType]->SetPoint(nPointsXY, hit.Position.X(), hit.Position.Y());
                    }
                }
            }
        }

        // Add hits to XZ projection (using RECONSTRUCTED hit positions)
        if (eventParticleGraphsXZ.find(particleType) != eventParticleGraphsXZ.end()) {
            for (int plane = 0; plane < 3; plane++) {
                for (UInt_t h = 0; h < particle->Hits[plane].size(); h++) {
                    AnaHitPD hit = particle->Hits[plane][h];
                    if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                        int nPointsXZ = eventParticleGraphsXZ[particleType]->GetN();
                        eventParticleGraphsXZ[particleType]->SetPoint(nPointsXZ, hit.Position.X(), hit.Position.Z());
                    }
                }
            }
        }
    }

    // Add vertices to graphs and mark parent particle positions
    for (const auto& vertex : vertices) {
        // Add vertex center
        int nPointsXY = vertexGraphXY->GetN();
        vertexGraphXY->SetPoint(nPointsXY, vertex.position[0], vertex.position[1]);

        int nPointsXZ = vertexGraphXZ->GetN();
        vertexGraphXZ->SetPoint(nPointsXZ, vertex.position[0], vertex.position[2]);

        // Only process the parent particle of this vertex
        if (vertex.parent) {
            // Add parent particle reconstructed end position (red squares)
            int nParentXY = parentEndGraphXY->GetN();
            parentEndGraphXY->SetPoint(nParentXY, vertex.parent->PositionEnd[0], vertex.parent->PositionEnd[1]);

            int nParentXZ = parentEndGraphXZ->GetN();
            parentEndGraphXZ->SetPoint(nParentXZ, vertex.parent->PositionEnd[0], vertex.parent->PositionEnd[2]);

            // Find and add the actual last hit of this parent particle (orange triangles)
            double lastHitX = -999, lastHitY = -999, lastHitZ = -999;
            double minResidualRange = 999999;
            bool foundLastHit = false;

            // Search through all planes to find the last hit of this parent using ResidualRange
            for (int plane = 0; plane < 3; plane++) {
                if (vertex.parent->Hits[plane].size() > 0) {
                    // Find the hit with the minimum ResidualRange (closest to end of track)
                    for (UInt_t h = 0; h < vertex.parent->Hits[plane].size(); h++) {
                        AnaHitPD hit = vertex.parent->Hits[plane][h];
                        if (hit.Position.X() > -900 && hit.Position.Y() > -900 && hit.ResidualRange > 0) {
                            if (!foundLastHit || hit.ResidualRange < minResidualRange) {
                                lastHitX = hit.Position.X();
                                lastHitY = hit.Position.Y();
                                lastHitZ = hit.Position.Z();
                                minResidualRange = hit.ResidualRange;
                                foundLastHit = true;
                            }
                        }
                    }
                }
            }

            // Add the last hit position if found
            if (foundLastHit) {
                int nLastHitXY = parentLastHitGraphXY->GetN();
                parentLastHitGraphXY->SetPoint(nLastHitXY, lastHitX, lastHitY);

                int nLastHitXZ = parentLastHitGraphXZ->GetN();
                parentLastHitGraphXZ->SetPoint(nLastHitXZ, lastHitX, lastHitZ);

                std::cout << "  Parent " << vertex.parent->GetTrueParticle()->PDG
                          << " (ID=" << vertex.parent->GetTrueParticle()->ID
                          << "): RecoEnd=(" << vertex.parent->PositionEnd[0]
                          << ", " << vertex.parent->PositionEnd[1]
                          << ", " << vertex.parent->PositionEnd[2]
                          << "), LastHit=(" << lastHitX
                          << ", " << lastHitY
                          << ", " << lastHitZ << "), ResRange=" << minResidualRange << std::endl;
            }
        }

        // Add first hits of daughter particles (all particles except the parent)
        for (size_t i = 1; i < vertex.particles.size(); i++) { // Skip index 0 (parent)
            AnaParticlePD* daughter = vertex.particles[i];
            AnaTrueParticleB* daughterTrue = daughter->GetTrueParticle();
            if (!daughterTrue) continue;

            // Find the first hit of this daughter particle
            double firstHitX = -999, firstHitY = -999, firstHitZ = -999;
            double maxResidualRange = -999;
            bool foundFirstHit = false;

            // Search through all planes to find the first hit using ResidualRange
            for (int plane = 0; plane < 3; plane++) {
                if (daughter->Hits[plane].size() > 0) {
                    // Find the hit with the maximum ResidualRange (closest to start of track)
                    for (UInt_t h = 0; h < daughter->Hits[plane].size(); h++) {
                        AnaHitPD hit = daughter->Hits[plane][h];
                        if (hit.Position.X() > -900 && hit.Position.Y() > -900 && hit.ResidualRange > 0) {
                            if (!foundFirstHit || hit.ResidualRange > maxResidualRange) {
                                firstHitX = hit.Position.X();
                                firstHitY = hit.Position.Y();
                                firstHitZ = hit.Position.Z();
                                maxResidualRange = hit.ResidualRange;
                                foundFirstHit = true;
                            }
                        }
                    }
                }
            }

            // Add the first hit position if found
            if (foundFirstHit) {
                int nFirstHitXY = daughterFirstHitGraphXY->GetN();
                daughterFirstHitGraphXY->SetPoint(nFirstHitXY, firstHitX, firstHitY);

                int nFirstHitXZ = daughterFirstHitGraphXZ->GetN();
                daughterFirstHitGraphXZ->SetPoint(nFirstHitXZ, firstHitX, firstHitZ);

                std::cout << "    Daughter " << daughterTrue->PDG
                          << " (ID=" << daughterTrue->ID
                          << "): FirstHit=(" << firstHitX
                          << ", " << firstHitY
                          << ", " << firstHitZ << "), ResRange=" << maxResidualRange << std::endl;
            }

            // Add the reconstructed start position of this daughter particle
            int nStartXY = daughterStartGraphXY->GetN();
            daughterStartGraphXY->SetPoint(nStartXY, daughter->PositionStart[0], daughter->PositionStart[1]);

            int nStartXZ = daughterStartGraphXZ->GetN();
            daughterStartGraphXZ->SetPoint(nStartXZ, daughter->PositionStart[0], daughter->PositionStart[2]);

            std::cout << "    Daughter " << daughterTrue->PDG
                      << " (ID=" << daughterTrue->ID
                      << "): RecoStart=(" << daughter->PositionStart[0]
                      << ", " << daughter->PositionStart[1]
                      << ", " << daughter->PositionStart[2] << ")" << std::endl;
        }
    }

    // Add particle ID labels to tracks
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue;

        // Find a representative point for each particle track (use first valid hit)
        double labelX = -999, labelY = -999, labelZ = -999;
        bool foundLabelPoint = false;

        for (int plane = 0; plane < 3 && !foundLabelPoint; plane++) {
            if (particle->Hits[plane].size() > 0) {
                for (UInt_t h = 0; h < particle->Hits[plane].size() && !foundLabelPoint; h++) {
                    AnaHitPD hit = particle->Hits[plane][h];
                    if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                        labelX = hit.Position.X();
                        labelY = hit.Position.Y();
                        labelZ = hit.Position.Z();
                        foundLabelPoint = true;
                    }
                }
            }
        }

        if (foundLabelPoint) {
            // Create text labels for XY projection
            std::string labelText = "ID:" + std::to_string(truePart->ID);
            TText* labelXY = new TText(labelX + 5, labelY + 5, labelText.c_str());
            labelXY->SetTextSize(0.02);
            labelXY->SetTextColor(kBlack);
            labelXY->SetTextFont(42);
            particleLabelsXY.push_back(labelXY);

            // Create text labels for XZ projection
            TText* labelXZ = new TText(labelX + 5, labelZ + 5, labelText.c_str());
            labelXZ->SetTextSize(0.02);
            labelXZ->SetTextColor(kBlack);
            labelXZ->SetTextFont(42);
            particleLabelsXZ.push_back(labelXZ);
        }
    }

    // Second loop: Find and draw K0 trajectories using the ParentIDs from first loop
    std::cout << "  Second loop: Finding K0 particles using " << k0ParentIDs.size() << " ParentIDs..." << std::endl;

    // Access all true particles in the spill
    std::vector<AnaTrueParticleB*>& allTrueParticles = spill->TrueParticles;
    Int_t nTrueParticles = spill->TrueParticles.size();

    std::cout << "  Total true particles in event: " << nTrueParticles << std::endl;

    // Find K0 particles by matching ParentIDs
    for (const auto& k0ParentID : k0ParentIDs) {
        for (const auto& truePart : allTrueParticles) {
            if (!truePart) continue;

            // Check if this is the K0 particle we're looking for
            if (truePart->ID == k0ParentID && truePart->PDG == 310) {
                std::cout << "  Found K0 particle: ID=" << truePart->ID << ", PDG=" << truePart->PDG << std::endl;
                std::cout << "  K0 trajectory: TrueStart=(" << truePart->Position[0]
                          << ", " << truePart->Position[1]
                          << ", " << truePart->Position[2]
                          << "), TrueEnd=(" << truePart->PositionEnd[0]
                          << ", " << truePart->PositionEnd[1]
                          << ", " << truePart->PositionEnd[2] << ")" << std::endl;

                // Create K0 trajectory lines using true start and end positions
                TLine* k0TrajectoryXY = new TLine(truePart->Position[0], truePart->Position[1],
                                                 truePart->PositionEnd[0], truePart->PositionEnd[1]);
                k0TrajectoryXY->SetLineColor(kViolet);
                k0TrajectoryXY->SetLineWidth(3);
                k0TrajectoryXY->SetLineStyle(1); // Solid line
                k0TrajectoriesXY.push_back(k0TrajectoryXY);

                TLine* k0TrajectoryXZ = new TLine(truePart->Position[0], truePart->Position[2],
                                                 truePart->PositionEnd[0], truePart->PositionEnd[2]);
                k0TrajectoryXZ->SetLineColor(kViolet);
                k0TrajectoryXZ->SetLineWidth(3);
                k0TrajectoryXZ->SetLineStyle(1); // Solid line
                k0TrajectoriesXZ.push_back(k0TrajectoryXZ);

                break; // Found this K0, move to next ParentID
            }
        }
    }

    // Draw XY projection (left subplot)
    c1->cd(1);
    std::string histTitleXY = "Event " + std::to_string(eventNumber) + " - XY Projection; X [cm]; Y [cm]";
    TH2F *dummyXY = new TH2F("dummyXY", histTitleXY.c_str(), 1000, -360, 360, 1000, 0, 700);
    dummyXY->SetStats(0); // Disable statistics box
    dummyXY->Draw();

    // Draw all particle graphs in XY projection
    for (auto& pair : eventParticleGraphsXY) {
        if (pair.second->GetN() > 0) {
            pair.second->Draw("P SAME");
        }
    }

    // Draw vertices in XY projection
    if (vertexGraphXY->GetN() > 0) {
        vertexGraphXY->Draw("P SAME");
    }

    // Draw parent particle end positions in XY projection
    if (parentEndGraphXY->GetN() > 0) {
        parentEndGraphXY->Draw("P SAME");
    }

    // Draw parent particle last hits in XY projection
    if (parentLastHitGraphXY->GetN() > 0) {
        parentLastHitGraphXY->Draw("P SAME");
    }

    // Draw daughter particle first hits in XY projection
    if (daughterFirstHitGraphXY->GetN() > 0) {
        daughterFirstHitGraphXY->Draw("P SAME");
    }

    // Draw daughter particle start positions in XY projection
    if (daughterStartGraphXY->GetN() > 0) {
        daughterStartGraphXY->Draw("P SAME");
    }

    // Draw vertex circles (30cm radius) in XY projection
    for (const auto& vertex : vertices) {
        TEllipse* vertexCircleXY = new TEllipse(vertex.position[0], vertex.position[1], 30.0, 30.0);
        vertexCircleXY->SetFillStyle(0); // Hollow circle
        vertexCircleXY->SetLineColor(kBlue);
        vertexCircleXY->SetLineWidth(2);
        vertexCircleXY->SetLineStyle(2); // Dashed line
        vertexCircleXY->Draw();
        vertexCirclesXY.push_back(vertexCircleXY);
    }

    // Draw particle ID labels in XY projection
    for (auto* label : particleLabelsXY) {
        if (label) label->Draw("SAME");
    }

    // Draw K0 trajectories in XY projection
    for (auto* k0Trajectory : k0TrajectoriesXY) {
        if (k0Trajectory) k0Trajectory->Draw("SAME");
    }


    // Draw beam particle trajectory in XY projection
    TLine* beamTrajectoryXY = new TLine(beamTrue->Position[0], beamTrue->Position[1],
                                       beamTrue->PositionEnd[0], beamTrue->PositionEnd[1]);
    beamTrajectoryXY->SetLineColor(kBlack);
    beamTrajectoryXY->SetLineWidth(3);
    beamTrajectoryXY->SetLineStyle(1);
    beamTrajectoryXY->Draw();

    // Draw XZ projection (right subplot)
    c1->cd(2);
    std::string histTitleXZ = "Event " + std::to_string(eventNumber) + " - XZ Projection; X [cm]; Z [cm]";
    TH2F *dummyXZ = new TH2F("dummyXZ", histTitleXZ.c_str(), 1000, -360, 360, 1000, 0, 700);
    dummyXZ->SetStats(0); // Disable statistics box
    dummyXZ->Draw();

    // Draw all particle graphs in XZ projection
    for (auto& pair : eventParticleGraphsXZ) {
        if (pair.second->GetN() > 0) {
            pair.second->Draw("P SAME");
        }
    }

    // Draw vertices in XZ projection
    if (vertexGraphXZ->GetN() > 0) {
        vertexGraphXZ->Draw("P SAME");
    }

    // Draw parent particle end positions in XZ projection
    if (parentEndGraphXZ->GetN() > 0) {
        parentEndGraphXZ->Draw("P SAME");
    }

    // Draw parent particle last hits in XZ projection
    if (parentLastHitGraphXZ->GetN() > 0) {
        parentLastHitGraphXZ->Draw("P SAME");
    }

    // Draw daughter particle first hits in XZ projection
    if (daughterFirstHitGraphXZ->GetN() > 0) {
        daughterFirstHitGraphXZ->Draw("P SAME");
    }

    // Draw daughter particle start positions in XZ projection
    if (daughterStartGraphXZ->GetN() > 0) {
        daughterStartGraphXZ->Draw("P SAME");
    }

    // Draw vertex circles (30cm radius) in XZ projection
    for (const auto& vertex : vertices) {
        TEllipse* vertexCircleXZ = new TEllipse(vertex.position[0], vertex.position[2], 30.0, 30.0);
        vertexCircleXZ->SetFillStyle(0); // Hollow circle
        vertexCircleXZ->SetLineColor(kBlue);
        vertexCircleXZ->SetLineWidth(2);
        vertexCircleXZ->SetLineStyle(2); // Dashed line
        vertexCircleXZ->Draw();
        vertexCirclesXZ.push_back(vertexCircleXZ);
    }

    // Draw particle ID labels in XZ projection
    for (auto* label : particleLabelsXZ) {
        if (label) label->Draw("SAME");
    }

    // Draw K0 trajectories in XZ projection
    for (auto* k0Trajectory : k0TrajectoriesXZ) {
        if (k0Trajectory) k0Trajectory->Draw("SAME");
    }


    // Draw beam particle trajectory in XZ projection
    TLine* beamTrajectoryXZ = new TLine(beamTrue->Position[0], beamTrue->Position[2],
                                       beamTrue->PositionEnd[0], beamTrue->PositionEnd[2]);
    beamTrajectoryXZ->SetLineColor(kBlack);
    beamTrajectoryXZ->SetLineWidth(3);
    beamTrajectoryXZ->SetLineStyle(1);
    beamTrajectoryXZ->Draw();

    // Create legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.95, 0.95);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(1);

    // Add entries for vertices
    legend->AddEntry(vertexGraphXY, "Vertex Centers", "p");

    // Add entries for parent particle end positions
    legend->AddEntry(parentEndGraphXY, "Parent End Positions", "p");

    // Add entries for parent particle last hits
    legend->AddEntry(parentLastHitGraphXY, "Parent Last Hits", "p");

    // Add entries for daughter particle first hits
    legend->AddEntry(daughterFirstHitGraphXY, "Daughter First Hits", "p");

    // Add entries for daughter particle start positions
    legend->AddEntry(daughterStartGraphXY, "Daughter Start Positions", "p");

    // Add entry for vertex circles
    TLine* vertexCircleLegend = new TLine(0, 0, 1, 1);
    vertexCircleLegend->SetLineColor(kBlue);
    vertexCircleLegend->SetLineWidth(2);
    vertexCircleLegend->SetLineStyle(2);
    legend->AddEntry(vertexCircleLegend, "Vertex Region (30cm)", "l");

    // Add entries for beam particle
    TLine* beamLegendLine = new TLine(0, 0, 1, 1);
    beamLegendLine->SetLineColor(kBlack);
    beamLegendLine->SetLineWidth(3);
    beamLegendLine->SetLineStyle(1);
    legend->AddEntry(beamLegendLine, "Beam Particle", "l");

    TLine* k0LegendLine = new TLine(0, 0, 1, 1);
    k0LegendLine->SetLineColor(kViolet);
    k0LegendLine->SetLineWidth(3);
    k0LegendLine->SetLineStyle(1);
    legend->AddEntry(k0LegendLine, "K0 Trajectory", "l");


    // Add entries for particle types that have hits (using TRUE PDG codes for identification)
    for (auto& pair : eventParticleGraphsXY) {
        if (pair.second && pair.second->GetN() > 0) {
            std::string particleName = "PDG " + std::to_string(pair.first);
            legend->AddEntry(pair.second, particleName.c_str(), "p");
        }
    }

    legend->Draw();

    c1->Update();
    std::cout << "Canvas is ready for interaction. Press any key in the canvas to continue..." << std::endl;
    c1->WaitPrimitive("k");

    // Clean up graphs (following PreliminaryK0Selection pattern)
    for (auto& pair : eventParticleGraphs) {
        if (pair.second) delete pair.second;
    }
    for (auto& pair : eventParticleGraphsXY) {
        if (pair.second) delete pair.second;
    }
    for (auto& pair : eventParticleGraphsXZ) {
        if (pair.second) delete pair.second;
    }
    if (vertexGraphXY) delete vertexGraphXY;
    if (vertexGraphXZ) delete vertexGraphXZ;
    if (parentEndGraphXY) delete parentEndGraphXY;
    if (parentEndGraphXZ) delete parentEndGraphXZ;
    if (parentLastHitGraphXY) delete parentLastHitGraphXY;
    if (parentLastHitGraphXZ) delete parentLastHitGraphXZ;
    if (daughterFirstHitGraphXY) delete daughterFirstHitGraphXY;
    if (daughterFirstHitGraphXZ) delete daughterFirstHitGraphXZ;
    if (daughterStartGraphXY) delete daughterStartGraphXY;
    if (daughterStartGraphXZ) delete daughterStartGraphXZ;

    // Clean up vertex circles
    for (auto* circle : vertexCirclesXY) {
        delete circle;
    }
    for (auto* circle : vertexCirclesXZ) {
        delete circle;
    }

    // Clean up particle labels
    for (auto* label : particleLabelsXY) {
        if (label) delete label;
    }
    for (auto* label : particleLabelsXZ) {
        if (label) delete label;
    }

    // Clean up K0 trajectory lines
    for (auto* k0Trajectory : k0TrajectoriesXY) {
        if (k0Trajectory) delete k0Trajectory;
    }
    for (auto* k0Trajectory : k0TrajectoriesXZ) {
        if (k0Trajectory) delete k0Trajectory;
    }


    // Note: TLine objects (beam trajectories) drawn on canvas are automatically cleaned up by ROOT
    // when the canvas is destroyed, so we don't need to manually delete them
    // Also, we don't close the canvas - let it remain open like PreliminaryK0Selection

    std::cout << "Event " << eventNumber << " visualization completed." << std::endl;
}
