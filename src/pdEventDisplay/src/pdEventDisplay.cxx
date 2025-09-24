#include "pdEventDisplay.hxx"
#include "pdAnalysisUtils.hxx"
#include "TDirectory.h"
#include <cstdlib>

// Initialize static member
int pdEventDisplay::_eventsDisplayed = 0;

//********************************************************************
pdEventDisplay::pdEventDisplay() {
//********************************************************************
    _CreateEventDisplay = false;
    _SaveToRootFile = false;
    _OutputDirectory = "/plots";
    _MaxEventsToDisplay = 100;
    _EventPercentage = 100.0;
    _RequiredParticlePDGs.clear();
    _VertexRadius = 30.0;
    _MinVertexDaughters = 2;
}

//********************************************************************
pdEventDisplay::~pdEventDisplay() {
//********************************************************************
    // Cleanup if needed
}

//********************************************************************
bool pdEventDisplay::Initialize(bool createEventDisplay, bool saveToRootFile, const std::string& outputDirectory, int maxEventsToDisplay, double eventPercentage, const std::vector<int>& requiredParticlePDGs, double vertexRadius, int minVertexDaughters){
//********************************************************************

    // Set parameters
    _CreateEventDisplay = createEventDisplay;
    _SaveToRootFile = saveToRootFile;
    _OutputDirectory = outputDirectory;
    _MaxEventsToDisplay = maxEventsToDisplay;
    _EventPercentage = eventPercentage;
    _RequiredParticlePDGs = requiredParticlePDGs;
    _VertexRadius = vertexRadius;
    _MinVertexDaughters = minVertexDaughters;

    // Initialize particle colors
    InitializeParticleColors();

    return true;
}


//********************************************************************
void pdEventDisplay::InitializeParticleColors(){
//********************************************************************

    // Initialize particle colors based on PDG codes - each particle type gets a unique color
    _particleColors[2212] = kRed+1;    // Proton (darker red - very distinct)
    _particleColors[-2212] = kRed+1;   // Anti-proton (same as proton)
    _particleColors[211] = kBlue;      // Pi+ (blue)
    _particleColors[-211] = kPink;     // Pi- (pink)
    _particleColors[111] = kCyan;      // Pi0 (cyan)
    _particleColors[321] = kGreen;     // K+ (green)
    _particleColors[-321] = kOrange;   // K- (orange)
    _particleColors[310] = kMagenta;   // K0 (magenta)
    _particleColors[130] = kMagenta;   // K0_L (same as K0)
    _particleColors[311] = kMagenta;   // K0_S (same as K0)
    _particleColors[13] = kYellow;     // Mu- (yellow)
    _particleColors[-13] = kYellow;    // Mu+ (same as mu-)
    _particleColors[11] = kTeal;       // e- (teal)
    _particleColors[-11] = kTeal;      // e+ (same as e-)
    _particleColors[22] = kGray+2;     // gamma (dark gray)
    _particleColors[2112] = kSpring;   // Neutron (spring green)
    _particleColors[-2112] = kSpring;  // Anti-neutron (same as neutron)
}

//********************************************************************
void pdEventDisplay::CreateEventDisplay(AnaEventB& event, int eventNumber){
//********************************************************************

    // Check if we've reached the maximum number of events to display
    if (_eventsDisplayed >= _MaxEventsToDisplay) {
        return;
    }

    // Check if this event should be saved based on percentage
    if (!ShouldSaveThisEvent(eventNumber)) {
        return;
    }

    _eventsDisplayed++;

    // Create the event display
    DrawEventProjections(event, eventNumber);
}

//********************************************************************
void pdEventDisplay::DrawEventProjections(AnaEventB& event, int eventNumber){
//********************************************************************


    // Get all particles in the event
    AnaParticleB** particles = event.Particles;
    int nParticles = event.nParticles;

    // Convert to vector of AnaParticlePD* for easier handling
    std::vector<AnaParticlePD*> allParticles;
    for (int i = 0; i < nParticles; i++) {
        AnaParticlePD* particle = static_cast<AnaParticlePD*>(particles[i]);
        if (particle) {
            allParticles.push_back(particle);
        }
    }

    // Get beam particle for vertex creation
    AnaTrueParticleB* beamTrue = nullptr;
    if (event.Beam) {
        AnaBeam* beam = static_cast<AnaBeam*>(event.Beam);
        if (beam && beam->BeamParticle) {
            beamTrue = beam->BeamParticle->GetTrueParticle();
        }
    }

    // Create vertices using the common utility function
    const double maxVertexRadius = _VertexRadius; // Use the same radius as configured
    std::vector<AnaVertexPD*> reconVertices = pdAnaUtils::CreateReconstructedVertices(event, maxVertexRadius);

    // Convert AnaVertexPD* to DisplayVertex for the event display
    std::vector<DisplayVertex> vertices;
    for (const auto& reconVertex : reconVertices) {
        if (!reconVertex || reconVertex->NParticles != 3) continue; // Skip if not 1 parent + 2 daughters

        DisplayVertex displayVertex;

        // Set vertex position
        displayVertex.position[0] = reconVertex->Position[0];
        displayVertex.position[1] = reconVertex->Position[1];
        displayVertex.position[2] = reconVertex->Position[2];

        // Set vertex properties
        displayVertex.parent = static_cast<AnaParticlePD*>(reconVertex->Particles[0]); // First particle is parent
        displayVertex.nParticles = 2; // exactly 2 daughters

        // Convert particles vector
        displayVertex.particles.clear();
        for (int i = 0; i < reconVertex->NParticles; i++) {
            displayVertex.particles.push_back(static_cast<AnaParticlePD*>(reconVertex->Particles[i]));
        }

        displayVertex.isTrueVertex = true; // For display purposes
        displayVertex.vertexType = 1; // True vertex

        vertices.push_back(displayVertex);
    }

    // Check if event contains signal vertices (K0 → π+ π- decays)
    // Use the same enhanced logic as EventContainsSignalVertices
    // Store signal flags in an array to avoid overwriting
    std::vector<bool> vertexSignalFlags(reconVertices.size(), false);
    int signalVertexCount = 0;

    for (size_t vtxIdx = 0; vtxIdx < reconVertices.size(); vtxIdx++) {
        const auto& reconVertex = reconVertices[vtxIdx];
        if (!reconVertex) continue;

        // Check if this vertex has exactly 2 daughters (K0 → pi+ pi-)
        if (reconVertex->NParticles >= 3) { // parent + 2 daughters
            std::vector<AnaTrueParticleB*> daughterTrueParts;

            // Skip the first particle (parent) and collect the daughters
            for (int i = 1; i < reconVertex->NParticles; i++) {
                AnaParticlePD* recoPart = reconVertex->Particles[i];
                if (!recoPart) continue;

                // Get the associated true particle
                AnaTrueParticleB* truePart = recoPart->GetTrueParticle();
                if (!truePart) continue;

                daughterTrueParts.push_back(truePart);
            }

            // Check if we have exactly 2 daughters
            if (daughterTrueParts.size() == 2) {
                AnaTrueParticleB* daughter1 = daughterTrueParts[0];
                AnaTrueParticleB* daughter2 = daughterTrueParts[1];

                // Check if both daughters have the same K0 parent (PDG=310 and same ParentID)
                if (daughter1->ParentPDG == 310 && daughter2->ParentPDG == 310) {
                    // Check if they have the same parent particle ID (same K0)
                    if (daughter1->ParentID == daughter2->ParentID && daughter1->ParentID != 0) {
                        // Now check if they are specifically π+ and π- (signal definition)
                        if ((daughter1->PDG == 211 && daughter2->PDG == -211) ||
                            (daughter1->PDG == -211 && daughter2->PDG == 211)) {

                            // Additional condition: Check if the K0 parent has a K+ parent (PDG=321)
                            // and that this K+ is the same particle that originated the vertex

                            // Find the K0 parent particle in the true particles list
                            AnaTrueParticleB* k0Parent = nullptr;
                            AnaTrueParticleB** allTrueParticles = event.TrueParticles;
                            int nTrueParticles = event.nTrueParticles;

                            for (int j = 0; j < nTrueParticles; j++) {
                                if (allTrueParticles[j] && allTrueParticles[j]->ID == daughter1->ParentID) {
                                    k0Parent = allTrueParticles[j];
                                    break;
                                }
                            }

                            if (k0Parent) {
                                // Check if the K0 parent has a K+ parent (PDG=321)
                                if (k0Parent->ParentPDG == 321) {
                                    // Find the K+ parent particle
                                    AnaTrueParticleB* kPlusParent = nullptr;
                                    for (int k = 0; k < nTrueParticles; k++) {
                                        if (allTrueParticles[k] && allTrueParticles[k]->ID == k0Parent->ParentID) {
                                            kPlusParent = allTrueParticles[k];
                                            break;
                                        }
                                    }

                                    if (kPlusParent) {
                                        // Check if this K+ is the same particle that originated the vertex
                                        // The vertex parent should be the reconstructed particle corresponding to this K+
                                        AnaParticlePD* vertexParent = static_cast<AnaParticlePD*>(reconVertex->Particles[0]);
                                        if (vertexParent) {
                                            AnaTrueParticleB* vertexParentTrue = vertexParent->GetTrueParticle();
                                            if (vertexParentTrue && vertexParentTrue->ID == kPlusParent->ID) {
                                                vertexSignalFlags[vtxIdx] = true; // Mark this vertex as signal
                                                signalVertexCount++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Check if any vertex is classified as signal
    bool hasSignalVertices = false;
    for (bool isSignal : vertexSignalFlags) {
        if (isSignal) {
            hasSignalVertices = true;
            break; // Found at least one signal vertex
        }
    }


    // Create canvas with two subplots (XY and XZ projections) using event-specific name
    std::string canvasName = "c" + std::to_string(eventNumber);
    std::string canvasTitle = "Event Display - Event " + std::to_string(eventNumber) + " (Vertices: " + std::to_string(vertices.size()) + ")";
    if (hasSignalVertices) {
        canvasTitle += " [K0 Signal]";
    }
    TCanvas* c1 = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1200, 600);
    c1->Divide(2, 1);

    SetupCanvas(c1);

    // We'll find K0s later after nTrueParticles and allTrueParticles are defined

    // Create graphs for different particle types
    std::map<int, TGraph*> eventParticleGraphsXY;
    std::map<int, TGraph*> eventParticleGraphsXZ;

    // Initialize graphs for each particle type
    for (auto& colorPair : _particleColors) {
        int pdg = colorPair.first;

        TGraph* graphXY = new TGraph();
        graphXY->SetMarkerStyle(20);
        graphXY->SetMarkerSize(0.5);
        graphXY->SetMarkerColor(colorPair.second);
        eventParticleGraphsXY[pdg] = graphXY;

        TGraph* graphXZ = new TGraph();
        graphXZ->SetMarkerStyle(20);
        graphXZ->SetMarkerSize(0.5);
        graphXZ->SetMarkerColor(colorPair.second);
        eventParticleGraphsXZ[pdg] = graphXZ;
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
    // Note: Colors will be set dynamically based on parent particle type
    TGraph* parentEndGraphXY = new TGraph();
    parentEndGraphXY->SetMarkerStyle(22);
    parentEndGraphXY->SetMarkerSize(1.8);
    parentEndGraphXY->SetMarkerColor(kRed); // Default color, will be overridden

    TGraph* parentEndGraphXZ = new TGraph();
    parentEndGraphXZ->SetMarkerStyle(22);
    parentEndGraphXZ->SetMarkerSize(1.8);
    parentEndGraphXZ->SetMarkerColor(kRed); // Default color, will be overridden

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

    // Storage for K0 start and end position markers
    std::vector<TMarker*> k0StartMarkersXY;
    std::vector<TMarker*> k0EndMarkersXY;
    std::vector<TMarker*> k0StartMarkersXZ;
    std::vector<TMarker*> k0EndMarkersXZ;

    // Storage for K0 parent trajectories and markers
    std::vector<TLine*> k0ParentTrajectoriesXY;
    std::vector<TLine*> k0ParentTrajectoriesXZ;
    std::vector<TMarker*> k0ParentStartMarkersXY;
    std::vector<TMarker*> k0ParentEndMarkersXY;
    std::vector<TMarker*> k0ParentStartMarkersXZ;
    std::vector<TMarker*> k0ParentEndMarkersXZ;

    // Storage for K0 daughter trajectories and markers
    std::vector<TLine*> k0DaughterTrajectoriesXY;
    std::vector<TLine*> k0DaughterTrajectoriesXZ;
    std::vector<TMarker*> k0DaughterStartMarkersXY;
    std::vector<TMarker*> k0DaughterEndMarkersXY;
    std::vector<TMarker*> k0DaughterStartMarkersXZ;
    std::vector<TMarker*> k0DaughterEndMarkersXZ;

    // Storage for K0 process labels
    std::vector<TText*> k0ProcessLabelsXY;
    std::vector<TText*> k0ProcessLabelsXZ;

    // First, process true particles (especially K0 particles that might not be reconstructed)
    AnaTrueParticleB** trueParticles = event.TrueParticles;
    int nTrueParticles = event.nTrueParticles;

    for (int i = 0; i < nTrueParticles; i++) {
        AnaTrueParticleB* truePart = trueParticles[i];
        if (!truePart) continue;

        // Only process K0 particles (PDG = 310) for now
        if (truePart->PDG != 310) continue;

        // Create graphs for K0 particles if they don't exist
        int particleType = truePart->PDG;
        if (eventParticleGraphsXY.find(particleType) == eventParticleGraphsXY.end()) {
            int color = GetParticleColor(particleType);

            TGraph* graphXY = new TGraph();
            graphXY->SetMarkerStyle(24); // Different marker for true particles
            graphXY->SetMarkerSize(1.0);
            graphXY->SetMarkerColor(color);
            eventParticleGraphsXY[particleType] = graphXY;

            TGraph* graphXZ = new TGraph();
            graphXZ->SetMarkerStyle(24); // Different marker for true particles
            graphXZ->SetMarkerSize(1.0);
            graphXZ->SetMarkerColor(color);
            eventParticleGraphsXZ[particleType] = graphXZ;
        }

        // Add K0 start and end positions to graphs
        if (truePart->Position[0] > -900 && truePart->Position[1] > -900) {
            int nPointsXY = eventParticleGraphsXY[particleType]->GetN();
            eventParticleGraphsXY[particleType]->SetPoint(nPointsXY, truePart->Position[0], truePart->Position[1]);
        }

        if (truePart->Position[0] > -900 && truePart->Position[2] > -900) {
            int nPointsXZ = eventParticleGraphsXZ[particleType]->GetN();
            eventParticleGraphsXZ[particleType]->SetPoint(nPointsXZ, truePart->Position[0], truePart->Position[2]);
        }

        // Add K0 end positions if different from start
        if (truePart->PositionEnd[0] > -900 && truePart->PositionEnd[1] > -900 &&
            (truePart->Position[0] != truePart->PositionEnd[0] || truePart->Position[1] != truePart->PositionEnd[1])) {
            int nPointsXY = eventParticleGraphsXY[particleType]->GetN();
            eventParticleGraphsXY[particleType]->SetPoint(nPointsXY, truePart->PositionEnd[0], truePart->PositionEnd[1]);
        }

        if (truePart->PositionEnd[0] > -900 && truePart->PositionEnd[2] > -900 &&
            (truePart->Position[0] != truePart->PositionEnd[0] || truePart->Position[2] != truePart->PositionEnd[2])) {
            int nPointsXZ = eventParticleGraphsXZ[particleType]->GetN();
            eventParticleGraphsXZ[particleType]->SetPoint(nPointsXZ, truePart->PositionEnd[0], truePart->PositionEnd[2]);
        }
    }

    // Fill XY and XZ graphs with hit data (using TRUE PDG for identification, RECO positions for display)
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue; // Skip particles without true information

        // Use TRUE particle PDG for identification
        int particleType = truePart->PDG;

        // Create graphs for this particle type if they don't exist
        if (eventParticleGraphsXY.find(particleType) == eventParticleGraphsXY.end()) {
            // Get color for this particle type (use default if not predefined)
            int color = GetParticleColor(particleType);

            TGraph* graphXY = new TGraph();
            graphXY->SetMarkerStyle(20);
            graphXY->SetMarkerSize(0.5);
            graphXY->SetMarkerColor(color);
            eventParticleGraphsXY[particleType] = graphXY;

            TGraph* graphXZ = new TGraph();
            graphXZ->SetMarkerStyle(20);
            graphXZ->SetMarkerSize(0.5);
            graphXZ->SetMarkerColor(color);
            eventParticleGraphsXZ[particleType] = graphXZ;
        }

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
                    if (hit.Position.X() > -900 && hit.Position.Z() > -900) {
                        int nPointsXZ = eventParticleGraphsXZ[particleType]->GetN();
                        eventParticleGraphsXZ[particleType]->SetPoint(nPointsXZ, hit.Position.X(), hit.Position.Z());
                    }
                }
            }
        }

        // Add particle ID labels
        if (eventParticleGraphsXY.find(particleType) != eventParticleGraphsXY.end() &&
            eventParticleGraphsXY[particleType]->GetN() > 0) {
            // Use the first hit position for labeling
            double labelX, labelY, labelZ;
            bool foundFirstHit = false;

            for (int plane = 0; plane < 3 && !foundFirstHit; plane++) {
                for (UInt_t h = 0; h < particle->Hits[plane].size() && !foundFirstHit; h++) {
                    AnaHitPD hit = particle->Hits[plane][h];
                    if (hit.Position.X() > -900 && hit.Position.Y() > -900 && hit.Position.Z() > -900) {
                        labelX = hit.Position.X();
                        labelY = hit.Position.Y();
                        labelZ = hit.Position.Z();
                        foundFirstHit = true;
                    }
                }
            }

            if (foundFirstHit) {
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
    }

    // Fill vertex information
    for (const auto& vertex : vertices) {
        // Add vertex center to graphs
        int nVertXY = vertexGraphXY->GetN();
        vertexGraphXY->SetPoint(nVertXY, vertex.position[0], vertex.position[1]);

        int nVertXZ = vertexGraphXZ->GetN();
        vertexGraphXZ->SetPoint(nVertXZ, vertex.position[0], vertex.position[2]);

        // Add parent particle end position if available
        if (vertex.parent) {
            // Get the color for the parent particle type
            AnaTrueParticleB* parentTrue = vertex.parent->GetTrueParticle();
            int parentColor = kRed; // Default color
            if (parentTrue && _particleColors.find(parentTrue->PDG) != _particleColors.end()) {
                parentColor = _particleColors[parentTrue->PDG];
            }

            int nParentEndXY = parentEndGraphXY->GetN();
            parentEndGraphXY->SetPoint(nParentEndXY, vertex.parent->PositionEnd[0], vertex.parent->PositionEnd[1]);

            int nParentEndXZ = parentEndGraphXZ->GetN();
            parentEndGraphXZ->SetPoint(nParentEndXZ, vertex.parent->PositionEnd[0], vertex.parent->PositionEnd[2]);

            // Set the color for the parent end position markers
            parentEndGraphXY->SetMarkerColor(parentColor);
            parentEndGraphXZ->SetMarkerColor(parentColor);

            // Add parent particle last hit (using ResidualRange to find the actual last hit)
            if (!vertex.parent->Hits[2].empty()) {
                AnaHitPD lastHit = vertex.parent->Hits[2][0]; // Initialize with first hit
                double minResidualRange = lastHit.ResidualRange;

                for (int plane = 0; plane < 3; plane++) {
                    for (UInt_t h = 0; h < vertex.parent->Hits[plane].size(); h++) {
                        AnaHitPD hit = vertex.parent->Hits[plane][h];
                        if (hit.ResidualRange < minResidualRange) {
                            minResidualRange = hit.ResidualRange;
                            lastHit = hit;
                        }
                    }
                }

                if (lastHit.Position.X() > -900 && lastHit.Position.Y() > -900) {
                    int nParentLastHitXY = parentLastHitGraphXY->GetN();
                    parentLastHitGraphXY->SetPoint(nParentLastHitXY, lastHit.Position.X(), lastHit.Position.Y());
                }

                if (lastHit.Position.X() > -900 && lastHit.Position.Z() > -900) {
                    int nParentLastHitXZ = parentLastHitGraphXZ->GetN();
                    parentLastHitGraphXZ->SetPoint(nParentLastHitXZ, lastHit.Position.X(), lastHit.Position.Z());
                }
            }
        }

        // Add daughter particle first hits and start positions
        for (const auto& daughter : vertex.particles) {
            if (daughter == vertex.parent) continue; // Skip the parent particle

            // Add daughter start position (reconstructed)
            int nDaughterStartXY = daughterStartGraphXY->GetN();
            daughterStartGraphXY->SetPoint(nDaughterStartXY, daughter->PositionStart[0], daughter->PositionStart[1]);

            int nDaughterStartXZ = daughterStartGraphXZ->GetN();
            daughterStartGraphXZ->SetPoint(nDaughterStartXZ, daughter->PositionStart[0], daughter->PositionStart[2]);

            // Add daughter first hit (using ResidualRange to find the actual first hit)
            if (!daughter->Hits[2].empty()) {
                AnaHitPD firstHit = daughter->Hits[2][0]; // Initialize with first hit
                double maxResidualRange = firstHit.ResidualRange;

                for (int plane = 0; plane < 3; plane++) {
                    for (UInt_t h = 0; h < daughter->Hits[plane].size(); h++) {
                        AnaHitPD hit = daughter->Hits[plane][h];
                        if (hit.ResidualRange > maxResidualRange) {
                            maxResidualRange = hit.ResidualRange;
                            firstHit = hit;
                        }
                    }
                }

                if (firstHit.Position.X() > -900 && firstHit.Position.Y() > -900) {
                    int nDaughterFirstHitXY = daughterFirstHitGraphXY->GetN();
                    daughterFirstHitGraphXY->SetPoint(nDaughterFirstHitXY, firstHit.Position.X(), firstHit.Position.Y());
                }

                if (firstHit.Position.X() > -900 && firstHit.Position.Z() > -900) {
                    int nDaughterFirstHitXZ = daughterFirstHitGraphXZ->GetN();
                    daughterFirstHitGraphXZ->SetPoint(nDaughterFirstHitXZ, firstHit.Position.X(), firstHit.Position.Z());
                }
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
        TEllipse* vertexCircleXY = new TEllipse(vertex.position[0], vertex.position[1], _VertexRadius, _VertexRadius);
        vertexCircleXY->SetFillStyle(0); // Hollow circle
        vertexCircleXY->SetLineColor(kBlue);
        vertexCircleXY->SetLineWidth(2);
        vertexCircleXY->Draw();
        vertexCirclesXY.push_back(vertexCircleXY);
    }

    // Draw particle ID labels in XY projection
    for (auto* label : particleLabelsXY) {
        if (label) label->Draw("SAME");
    }

    // Second loop: Find and draw K0 trajectories using the ParentIDs from first loop
    // Access all true particles in the spill
    AnaTrueParticleB** allTrueParticles = event.TrueParticles;
    // Int_t nTrueParticles = event.nTrueParticles;

    // Draw ALL K0s in the event, not just those that are parents of pions in vertices
    std::vector<Int_t> k0ParentIDs;
    for (int i = 0; i < nTrueParticles; i++) {
        AnaTrueParticleB* truePart = allTrueParticles[i];
        if (!truePart) continue;

        // If this is a K0, add it to our list
        if (truePart->PDG == 310) {
            bool alreadyStored = false;
            for (const auto& existingID : k0ParentIDs) {
                if (existingID == truePart->ID) {
                    alreadyStored = true;
                    break;
                }
            }
            if (!alreadyStored) {
                k0ParentIDs.push_back(truePart->ID);
            }
        }
    }

    // std::cout << "Event " << eventNumber << ": Found " << k0ParentIDs.size() << " K0 particles to draw" << std::endl;

    // Find K0 particles by matching ParentIDs
    for (const auto& k0ParentID : k0ParentIDs) {
        for (int i = 0; i < nTrueParticles; i++) {
            AnaTrueParticleB* truePart = allTrueParticles[i];
            if (!truePart) continue;

            // Check if this is the K0 particle we're looking for
            if (truePart->ID == k0ParentID && truePart->PDG == 310) {
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

                // Create K0 start position markers
                TMarker* k0StartMarkerXY = new TMarker(truePart->Position[0], truePart->Position[1], kFullCircle);
                k0StartMarkerXY->SetMarkerColor(kViolet);
                k0StartMarkerXY->SetMarkerSize(1.0);
                k0StartMarkersXY.push_back(k0StartMarkerXY);

                TMarker* k0StartMarkerXZ = new TMarker(truePart->Position[0], truePart->Position[2], kFullCircle);
                k0StartMarkerXZ->SetMarkerColor(kViolet);
                k0StartMarkerXZ->SetMarkerSize(1.0);
                k0StartMarkersXZ.push_back(k0StartMarkerXZ);

                // Create K0 end position markers
                TMarker* k0EndMarkerXY = new TMarker(truePart->PositionEnd[0], truePart->PositionEnd[1], kFullSquare);
                k0EndMarkerXY->SetMarkerColor(kViolet);
                k0EndMarkerXY->SetMarkerSize(1.0);
                k0EndMarkersXY.push_back(k0EndMarkerXY);

                TMarker* k0EndMarkerXZ = new TMarker(truePart->PositionEnd[0], truePart->PositionEnd[2], kFullSquare);
                k0EndMarkerXZ->SetMarkerColor(kViolet);
                k0EndMarkerXZ->SetMarkerSize(1.0);
                k0EndMarkersXZ.push_back(k0EndMarkerXZ);

                // Cast to AnaTrueParticlePD to access the Daughters field
                AnaTrueParticlePD* truePartPD = static_cast<AnaTrueParticlePD*>(truePart);
                int nDaughters = 0;

                if (truePartPD && !truePartPD->Daughters.empty()) {
                    // Use the Daughters field to get the daughter IDs
                    nDaughters = truePartPD->Daughters.size();

                    // Now find and draw the actual daughter particles using the Daughters field
                    for (size_t d = 0; d < truePartPD->Daughters.size(); d++) {
                        int daughterID = truePartPD->Daughters[d];

                        // Find the daughter particle by ID
                        AnaTrueParticleB* daughterPart = nullptr;
                        for (int j = 0; j < nTrueParticles; j++) {
                            if (allTrueParticles[j] && allTrueParticles[j]->ID == daughterID) {
                                daughterPart = allTrueParticles[j];
                                break;
                            }
                        }

                        if (!daughterPart) continue; // Daughter not found in event

                        // Get color for this daughter particle type
                        int daughterColor = GetParticleColor(daughterPart->PDG);

                        // Check if this daughter has been reconstructed
                        bool isReconstructed = false;
                        for (const auto& particle : allParticles) {
                            AnaTrueParticleB* truePartReco = particle->GetTrueParticle();
                            if (truePartReco && truePartReco->ID == daughterPart->ID) {
                                isReconstructed = true;
                                break;
                            }
                        }

                        // Create daughter trajectory lines using true start and end positions
                        TLine* daughterTrajectoryXY = new TLine(daughterPart->Position[0], daughterPart->Position[1],
                                                               daughterPart->PositionEnd[0], daughterPart->PositionEnd[1]);
                        daughterTrajectoryXY->SetLineColor(daughterColor);
                        daughterTrajectoryXY->SetLineWidth(2);
                        // Use solid line for reconstructed daughters, dashed for non-reconstructed
                        daughterTrajectoryXY->SetLineStyle(isReconstructed ? 1 : 2);
                        k0DaughterTrajectoriesXY.push_back(daughterTrajectoryXY);

                        TLine* daughterTrajectoryXZ = new TLine(daughterPart->Position[0], daughterPart->Position[2],
                                                               daughterPart->PositionEnd[0], daughterPart->PositionEnd[2]);
                        daughterTrajectoryXZ->SetLineColor(daughterColor);
                        daughterTrajectoryXZ->SetLineWidth(2);
                        // Use solid line for reconstructed daughters, dashed for non-reconstructed
                        daughterTrajectoryXZ->SetLineStyle(isReconstructed ? 1 : 2);
                        k0DaughterTrajectoriesXZ.push_back(daughterTrajectoryXZ);

                        // Create daughter start position markers
                        TMarker* daughterStartMarkerXY = new TMarker(daughterPart->Position[0], daughterPart->Position[1], kFullCircle);
                        daughterStartMarkerXY->SetMarkerColor(daughterColor);
                        daughterStartMarkerXY->SetMarkerSize(0.8);
                        k0DaughterStartMarkersXY.push_back(daughterStartMarkerXY);

                        TMarker* daughterStartMarkerXZ = new TMarker(daughterPart->Position[0], daughterPart->Position[2], kFullCircle);
                        daughterStartMarkerXZ->SetMarkerColor(daughterColor);
                        daughterStartMarkerXZ->SetMarkerSize(0.8);
                        k0DaughterStartMarkersXZ.push_back(daughterStartMarkerXZ);

                        // Create daughter end position markers
                        TMarker* daughterEndMarkerXY = new TMarker(daughterPart->PositionEnd[0], daughterPart->PositionEnd[1], kFullSquare);
                        daughterEndMarkerXY->SetMarkerColor(daughterColor);
                        daughterEndMarkerXY->SetMarkerSize(0.8);
                        k0DaughterEndMarkersXY.push_back(daughterEndMarkerXY);

                        TMarker* daughterEndMarkerXZ = new TMarker(daughterPart->PositionEnd[0], daughterPart->PositionEnd[2], kFullSquare);
                        daughterEndMarkerXZ->SetMarkerColor(daughterColor);
                        daughterEndMarkerXZ->SetMarkerSize(0.8);
                        k0DaughterEndMarkersXZ.push_back(daughterEndMarkerXZ);
                    }
                } else {
                    // Fallback to the old method if Daughters field is not available
                    for (int j = 0; j < nTrueParticles; j++) {
                        AnaTrueParticleB* daughterPart = allTrueParticles[j];
                        if (!daughterPart) continue;
                        if (daughterPart->ParentID == truePart->ID) {
                            nDaughters++;

                            // Get color for this daughter particle type
                            int daughterColor = GetParticleColor(daughterPart->PDG);

                            // Check if this daughter has been reconstructed
                            bool isReconstructed = false;
                            for (const auto& particle : allParticles) {
                                AnaTrueParticleB* truePartReco = particle->GetTrueParticle();
                                if (truePartReco && truePartReco->ID == daughterPart->ID) {
                                    isReconstructed = true;
                                    break;
                                }
                            }

                            // Create daughter trajectory lines using true start and end positions
                            TLine* daughterTrajectoryXY = new TLine(daughterPart->Position[0], daughterPart->Position[1],
                                                                   daughterPart->PositionEnd[0], daughterPart->PositionEnd[1]);
                            daughterTrajectoryXY->SetLineColor(daughterColor);
                            daughterTrajectoryXY->SetLineWidth(2);
                            // Use solid line for reconstructed daughters, dashed for non-reconstructed
                            daughterTrajectoryXY->SetLineStyle(isReconstructed ? 1 : 2);
                            k0DaughterTrajectoriesXY.push_back(daughterTrajectoryXY);

                            TLine* daughterTrajectoryXZ = new TLine(daughterPart->Position[0], daughterPart->Position[2],
                                                                   daughterPart->PositionEnd[0], daughterPart->PositionEnd[2]);
                            daughterTrajectoryXZ->SetLineColor(daughterColor);
                            daughterTrajectoryXZ->SetLineWidth(2);
                            // Use solid line for reconstructed daughters, dashed for non-reconstructed
                            daughterTrajectoryXZ->SetLineStyle(isReconstructed ? 1 : 2);
                            k0DaughterTrajectoriesXZ.push_back(daughterTrajectoryXZ);

                            // Create daughter start position markers
                            TMarker* daughterStartMarkerXY = new TMarker(daughterPart->Position[0], daughterPart->Position[1], kFullCircle);
                            daughterStartMarkerXY->SetMarkerColor(daughterColor);
                            daughterStartMarkerXY->SetMarkerSize(0.8);
                            k0DaughterStartMarkersXY.push_back(daughterStartMarkerXY);

                            TMarker* daughterStartMarkerXZ = new TMarker(daughterPart->Position[0], daughterPart->Position[2], kFullCircle);
                            daughterStartMarkerXZ->SetMarkerColor(daughterColor);
                            daughterStartMarkerXZ->SetMarkerSize(0.8);
                            k0DaughterStartMarkersXZ.push_back(daughterStartMarkerXZ);

                            // Create daughter end position markers
                            TMarker* daughterEndMarkerXY = new TMarker(daughterPart->PositionEnd[0], daughterPart->PositionEnd[1], kFullSquare);
                            daughterEndMarkerXY->SetMarkerColor(daughterColor);
                            daughterEndMarkerXY->SetMarkerSize(0.8);
                            k0DaughterEndMarkersXY.push_back(daughterEndMarkerXY);

                            TMarker* daughterEndMarkerXZ = new TMarker(daughterPart->PositionEnd[0], daughterPart->PositionEnd[2], kFullSquare);
                            daughterEndMarkerXZ->SetMarkerColor(daughterColor);
                            daughterEndMarkerXZ->SetMarkerSize(0.8);
                            k0DaughterEndMarkersXZ.push_back(daughterEndMarkerXZ);
                        }
                    }
                }

                // Create process label next to the K0 with daughter count
                std::string processName = ConvertProcessToString(truePart->ProcessEnd);
                std::string processLabel = processName + ": " + std::to_string(nDaughters);
                TText* k0ProcessLabelXY = new TText(truePart->PositionEnd[0] + 5, truePart->PositionEnd[1] + 5, processLabel.c_str());
                k0ProcessLabelXY->SetTextColor(kViolet);
                k0ProcessLabelXY->SetTextSize(0.02);
                k0ProcessLabelsXY.push_back(k0ProcessLabelXY);

                TText* k0ProcessLabelXZ = new TText(truePart->PositionEnd[0] + 5, truePart->PositionEnd[2] + 5, processLabel.c_str());
                k0ProcessLabelXZ->SetTextColor(kViolet);
                k0ProcessLabelXZ->SetTextSize(0.02);
                k0ProcessLabelsXZ.push_back(k0ProcessLabelXZ);

                // Draw the true parent of the K0 if it exists
                if (truePart->ParentID != 0) {
                    // Find the parent particle by matching ParentID
                    for (int k = 0; k < nTrueParticles; k++) {
                        AnaTrueParticleB* parentPart = allTrueParticles[k];
                        if (!parentPart) continue;

                        // Check if this is the parent of the current K0
                        if (parentPart->ID == truePart->ParentID) {
                            // Get color for the parent particle type
                            int parentColor = GetParticleColor(parentPart->PDG);

                            // Create parent trajectory lines using true start and end positions
                            TLine* parentTrajectoryXY = new TLine(parentPart->Position[0], parentPart->Position[1],
                                                                 parentPart->PositionEnd[0], parentPart->PositionEnd[1]);
                            parentTrajectoryXY->SetLineColor(parentColor);
                            parentTrajectoryXY->SetLineWidth(2);
                            parentTrajectoryXY->SetLineStyle(1); // Solid line
                            k0ParentTrajectoriesXY.push_back(parentTrajectoryXY);

                            TLine* parentTrajectoryXZ = new TLine(parentPart->Position[0], parentPart->Position[2],
                                                                 parentPart->PositionEnd[0], parentPart->PositionEnd[2]);
                            parentTrajectoryXZ->SetLineColor(parentColor);
                            parentTrajectoryXZ->SetLineWidth(2);
                            parentTrajectoryXZ->SetLineStyle(1); // Solid line
                            k0ParentTrajectoriesXZ.push_back(parentTrajectoryXZ);

                            // Create parent start position markers
                            TMarker* parentStartMarkerXY = new TMarker(parentPart->Position[0], parentPart->Position[1], kFullCircle);
                            parentStartMarkerXY->SetMarkerColor(parentColor);
                            parentStartMarkerXY->SetMarkerSize(0.8);
                            k0ParentStartMarkersXY.push_back(parentStartMarkerXY);

                            TMarker* parentStartMarkerXZ = new TMarker(parentPart->Position[0], parentPart->Position[2], kFullCircle);
                            parentStartMarkerXZ->SetMarkerColor(parentColor);
                            parentStartMarkerXZ->SetMarkerSize(0.8);
                            k0ParentStartMarkersXZ.push_back(parentStartMarkerXZ);

                            // Create parent end position markers
                            TMarker* parentEndMarkerXY = new TMarker(parentPart->PositionEnd[0], parentPart->PositionEnd[1], kFullSquare);
                            parentEndMarkerXY->SetMarkerColor(parentColor);
                            parentEndMarkerXY->SetMarkerSize(0.8);
                            k0ParentEndMarkersXY.push_back(parentEndMarkerXY);

                            TMarker* parentEndMarkerXZ = new TMarker(parentPart->PositionEnd[0], parentPart->PositionEnd[2], kFullSquare);
                            parentEndMarkerXZ->SetMarkerColor(parentColor);
                            parentEndMarkerXZ->SetMarkerSize(0.8);
                            k0ParentEndMarkersXZ.push_back(parentEndMarkerXZ);

                            break; // Found the parent, no need to continue searching
                        }
                    }
                }


                break; // Found this K0, move to next ParentID
            }
        }
    }

    // Draw K0 trajectories in XY projection
    for (auto* k0Trajectory : k0TrajectoriesXY) {
        if (k0Trajectory) k0Trajectory->Draw("SAME");
    }

    // Draw K0 parent trajectories in XY projection
    for (auto* k0ParentTrajectory : k0ParentTrajectoriesXY) {
        if (k0ParentTrajectory) k0ParentTrajectory->Draw("SAME");
    }

    // Draw K0 start markers in XY projection
    for (auto* k0StartMarker : k0StartMarkersXY) {
        if (k0StartMarker) k0StartMarker->Draw("SAME");
    }

    // Draw K0 parent start markers in XY projection
    for (auto* k0ParentStartMarker : k0ParentStartMarkersXY) {
        if (k0ParentStartMarker) k0ParentStartMarker->Draw("SAME");
    }

    // Draw K0 end markers in XY projection
    for (auto* k0EndMarker : k0EndMarkersXY) {
        if (k0EndMarker) k0EndMarker->Draw("SAME");
    }

    // Draw K0 parent end markers in XY projection
    for (auto* k0ParentEndMarker : k0ParentEndMarkersXY) {
        if (k0ParentEndMarker) k0ParentEndMarker->Draw("SAME");
    }

    // Draw K0 process labels in XY projection
    for (auto* k0ProcessLabel : k0ProcessLabelsXY) {
        if (k0ProcessLabel) k0ProcessLabel->Draw("SAME");
    }

    // Draw K0 daughter trajectories in XY projection
    for (auto* daughterTrajectory : k0DaughterTrajectoriesXY) {
        if (daughterTrajectory) daughterTrajectory->Draw("SAME");
    }

    // Draw K0 daughter start markers in XY projection
    for (auto* daughterStartMarker : k0DaughterStartMarkersXY) {
        if (daughterStartMarker) daughterStartMarker->Draw("SAME");
    }

    // Draw K0 daughter end markers in XY projection
    for (auto* daughterEndMarker : k0DaughterEndMarkersXY) {
        if (daughterEndMarker) daughterEndMarker->Draw("SAME");
    }

    // Draw XZ projection (right subplot)
    c1->cd(2);
    std::string histTitleXZ = "Event " + std::to_string(eventNumber) + " - XZ Projection; X [cm]; Z [cm]";
    TH2F *dummyXZ = new TH2F("dummyXZ", histTitleXZ.c_str(), 1000, -360, 360, 1000, -600, 600);
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
        TEllipse* vertexCircleXZ = new TEllipse(vertex.position[0], vertex.position[2], _VertexRadius, _VertexRadius);
        vertexCircleXZ->SetFillStyle(0); // Hollow circle
        vertexCircleXZ->SetLineColor(kBlue);
        vertexCircleXZ->SetLineWidth(2);
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

    // Draw K0 parent trajectories in XZ projection
    for (auto* k0ParentTrajectory : k0ParentTrajectoriesXZ) {
        if (k0ParentTrajectory) k0ParentTrajectory->Draw("SAME");
    }

    // Draw K0 start markers in XZ projection
    for (auto* k0StartMarker : k0StartMarkersXZ) {
        if (k0StartMarker) k0StartMarker->Draw("SAME");
    }

    // Draw K0 parent start markers in XZ projection
    for (auto* k0ParentStartMarker : k0ParentStartMarkersXZ) {
        if (k0ParentStartMarker) k0ParentStartMarker->Draw("SAME");
    }

    // Draw K0 end markers in XZ projection
    for (auto* k0EndMarker : k0EndMarkersXZ) {
        if (k0EndMarker) k0EndMarker->Draw("SAME");
    }

    // Draw K0 parent end markers in XZ projection
    for (auto* k0ParentEndMarker : k0ParentEndMarkersXZ) {
        if (k0ParentEndMarker) k0ParentEndMarker->Draw("SAME");
    }

    // Draw K0 process labels in XZ projection
    for (auto* k0ProcessLabel : k0ProcessLabelsXZ) {
        if (k0ProcessLabel) k0ProcessLabel->Draw("SAME");
    }

    // Draw K0 daughter trajectories in XZ projection
    for (auto* daughterTrajectory : k0DaughterTrajectoriesXZ) {
        if (daughterTrajectory) daughterTrajectory->Draw("SAME");
    }

    // Draw K0 daughter start markers in XZ projection
    for (auto* daughterStartMarker : k0DaughterStartMarkersXZ) {
        if (daughterStartMarker) daughterStartMarker->Draw("SAME");
    }

    // Draw K0 daughter end markers in XZ projection
    for (auto* daughterEndMarker : k0DaughterEndMarkersXZ) {
        if (daughterEndMarker) daughterEndMarker->Draw("SAME");
    }

    // Create legend
    TLegend* legend = new TLegend(0.65, 0.65, 0.95, 0.95);
    legend->SetNColumns(1);
    legend->SetTextSize(0.025);

    // Add legend entries for particles that have hits
    for (auto& pair : eventParticleGraphsXY) {
        if (pair.second->GetN() > 0) {
            std::string particleName = GetParticleName(pair.first);
            legend->AddEntry(pair.second, particleName.c_str(), "P");
        }
    }

    // Add special legend entries
    if (vertexGraphXY->GetN() > 0) {
        legend->AddEntry(vertexGraphXY, "Vertex Center", "P");
    }
    if (parentEndGraphXY->GetN() > 0) {
        legend->AddEntry(parentEndGraphXY, "Parent End Pos", "P");
    }
    if (parentLastHitGraphXY->GetN() > 0) {
        legend->AddEntry(parentLastHitGraphXY, "Parent Last Hit", "P");
    }
    if (daughterFirstHitGraphXY->GetN() > 0) {
        legend->AddEntry(daughterFirstHitGraphXY, "Daughter First Hit", "P");
    }
    if (daughterStartGraphXY->GetN() > 0) {
        legend->AddEntry(daughterStartGraphXY, "Daughter Start Pos", "P");
    }
    if (k0TrajectoriesXY.size() > 0) {
        legend->AddEntry(k0TrajectoriesXY[0], "K0 Trajectory", "L");
    }
    if (k0ParentTrajectoriesXY.size() > 0) {
        legend->AddEntry(k0ParentTrajectoriesXY[0], "K0 Parent Trajectory", "L");
    }
    if (k0StartMarkersXY.size() > 0) {
        legend->AddEntry(k0StartMarkersXY[0], "K0 Start", "P");
    }
    if (k0ParentStartMarkersXY.size() > 0) {
        legend->AddEntry(k0ParentStartMarkersXY[0], "K0 Parent Start", "P");
    }
    if (k0EndMarkersXY.size() > 0) {
        legend->AddEntry(k0EndMarkersXY[0], "K0 End", "P");
    }
    if (k0ParentEndMarkersXY.size() > 0) {
        legend->AddEntry(k0ParentEndMarkersXY[0], "K0 Parent End", "P");
    }

    // Add legend entries for K0 daughters (reconstructed and non-reconstructed)
    if (k0DaughterTrajectoriesXY.size() > 0) {
        // Check if we have both reconstructed and non-reconstructed daughters
        bool hasReconstructed = false;
        bool hasNonReconstructed = false;

        for (auto* trajectory : k0DaughterTrajectoriesXY) {
            if (trajectory->GetLineStyle() == 1) hasReconstructed = true;
            if (trajectory->GetLineStyle() == 2) hasNonReconstructed = true;
        }

        if (hasReconstructed && hasNonReconstructed) {
            // Create dummy lines for legend
            TLine* dummyReconstructed = new TLine();
            dummyReconstructed->SetLineStyle(1);
            dummyReconstructed->SetLineColor(kBlack);
            legend->AddEntry(dummyReconstructed, "K0 Daughters (Reconstructed)", "L");

            TLine* dummyNonReconstructed = new TLine();
            dummyNonReconstructed->SetLineStyle(2);
            dummyNonReconstructed->SetLineColor(kBlack);
            legend->AddEntry(dummyNonReconstructed, "K0 Daughters (Not Reconstructed)", "L");
        } else if (hasReconstructed) {
            TLine* dummyReconstructed = new TLine();
            dummyReconstructed->SetLineStyle(1);
            dummyReconstructed->SetLineColor(kBlack);
            legend->AddEntry(dummyReconstructed, "K0 Daughters (Reconstructed)", "L");
        } else if (hasNonReconstructed) {
            TLine* dummyNonReconstructed = new TLine();
            dummyNonReconstructed->SetLineStyle(2);
            dummyNonReconstructed->SetLineColor(kBlack);
            legend->AddEntry(dummyNonReconstructed, "K0 Daughters (Not Reconstructed)", "L");
        }
    }

    DrawLegend(legend);

    // Update canvas
    c1->Update();

    // Save canvas to ROOT file if enabled
    if (_SaveToRootFile) {
        SaveCanvasToRootFile(c1, eventNumber);
    }

    // Clean up graphs
    for (auto& pair : eventParticleGraphsXY) {
        if (pair.second) delete pair.second;
    }
    for (auto& pair : eventParticleGraphsXZ) {
        if (pair.second) delete pair.second;
    }

    // Clean up special graphs
    delete vertexGraphXY;
    delete vertexGraphXZ;
    delete parentEndGraphXY;
    delete parentEndGraphXZ;
    delete parentLastHitGraphXY;
    delete parentLastHitGraphXZ;
    delete daughterFirstHitGraphXY;
    delete daughterFirstHitGraphXZ;
    delete daughterStartGraphXY;
    delete daughterStartGraphXZ;

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

    // Clean up K0 parent trajectory lines
    for (auto* k0ParentTrajectory : k0ParentTrajectoriesXY) {
        if (k0ParentTrajectory) delete k0ParentTrajectory;
    }
    for (auto* k0ParentTrajectory : k0ParentTrajectoriesXZ) {
        if (k0ParentTrajectory) delete k0ParentTrajectory;
    }

    // Clean up K0 markers
    for (auto* k0StartMarker : k0StartMarkersXY) {
        if (k0StartMarker) delete k0StartMarker;
    }
    for (auto* k0EndMarker : k0EndMarkersXY) {
        if (k0EndMarker) delete k0EndMarker;
    }
    for (auto* k0StartMarker : k0StartMarkersXZ) {
        if (k0StartMarker) delete k0StartMarker;
    }
    for (auto* k0EndMarker : k0EndMarkersXZ) {
        if (k0EndMarker) delete k0EndMarker;
    }

    // Clean up K0 parent markers
    for (auto* k0ParentStartMarker : k0ParentStartMarkersXY) {
        if (k0ParentStartMarker) delete k0ParentStartMarker;
    }
    for (auto* k0ParentEndMarker : k0ParentEndMarkersXY) {
        if (k0ParentEndMarker) delete k0ParentEndMarker;
    }
    for (auto* k0ParentStartMarker : k0ParentStartMarkersXZ) {
        if (k0ParentStartMarker) delete k0ParentStartMarker;
    }
    for (auto* k0ParentEndMarker : k0ParentEndMarkersXZ) {
        if (k0ParentEndMarker) delete k0ParentEndMarker;
    }

    // Clean up K0 daughter trajectories
    for (auto* daughterTrajectory : k0DaughterTrajectoriesXY) {
        if (daughterTrajectory) delete daughterTrajectory;
    }
    for (auto* daughterTrajectory : k0DaughterTrajectoriesXZ) {
        if (daughterTrajectory) delete daughterTrajectory;
    }

    // Clean up K0 daughter markers
    for (auto* daughterStartMarker : k0DaughterStartMarkersXY) {
        if (daughterStartMarker) delete daughterStartMarker;
    }
    for (auto* daughterEndMarker : k0DaughterEndMarkersXY) {
        if (daughterEndMarker) delete daughterEndMarker;
    }
    for (auto* daughterStartMarker : k0DaughterStartMarkersXZ) {
        if (daughterStartMarker) delete daughterStartMarker;
    }
    for (auto* daughterEndMarker : k0DaughterEndMarkersXZ) {
        if (daughterEndMarker) delete daughterEndMarker;
    }

    // Clean up K0 process labels
    for (auto* k0ProcessLabel : k0ProcessLabelsXY) {
        if (k0ProcessLabel) delete k0ProcessLabel;
    }
    for (auto* k0ProcessLabel : k0ProcessLabelsXZ) {
        if (k0ProcessLabel) delete k0ProcessLabel;
    }

    // Note: We don't delete the canvas here as it might be displayed on screen
    // The canvas will be cleaned up by ROOT's garbage collection
}

//********************************************************************
void pdEventDisplay::SaveCanvasToRootFile(TCanvas* canvas, int eventNumber){
//********************************************************************

    // Set a more descriptive title that will show up in ROOT
    std::string title = canvas->GetTitle();
    canvas->SetTitle(title.c_str());

    // The canvas already has the correct name (c + eventNumber), so just write it directly
    canvas->Write();
}

//********************************************************************
bool pdEventDisplay::ShouldSaveThisEvent(int eventNumber) const {
//********************************************************************

    // If percentage is 100%, save all events
    if (_EventPercentage >= 100.0) {
        return true;
    }

    // If percentage is 0%, save no events
    if (_EventPercentage <= 0.0) {
        return false;
    }

    // Batch-based approach: save X events from every 100 events
    // Calculate which batch of 100 this event belongs to
    int batchNumber = eventNumber / 100;
    int eventInBatch = eventNumber % 100;

    // Calculate how many events to save from each batch
    int eventsToSavePerBatch = (int)(_EventPercentage / 100.0 * 100.0);

    // Use a deterministic selection within each batch based on event number
    // This ensures the same events are always selected for the same percentage
    int selectionSeed = batchNumber * 1000 + eventInBatch; // Create a unique seed for this event
    double eventRatio = (double)(selectionSeed % 100) / 100.0;
    double threshold = (double)eventsToSavePerBatch / 100.0;

    return eventRatio < threshold;
}

//********************************************************************
int pdEventDisplay::GetParticleColor(int pdg){
//********************************************************************

    auto it = _particleColors.find(pdg);
    if (it != _particleColors.end()) {
        return it->second;
    }
    return kBlack; // Default color
}

//********************************************************************
std::string pdEventDisplay::GetParticleName(int pdg){
//********************************************************************

    switch (pdg) {
        case 2212:  return "Proton";
        case -2212: return "Anti-proton";
        case 211:   return "Pi+";
        case -211:  return "Pi-";
        case 111:   return "Pi0";
        case 321:   return "K+";
        case -321:  return "K-";
        case 310:   return "K0";
        case 130:   return "K0_L";
        case 311:   return "K0_S";
        case 13:    return "Mu-";
        case -13:   return "Mu+";
        case 11:    return "e-";
        case -11:   return "e+";
        case 22:    return "Gamma";
        case 2112:  return "Neutron";
        case -2112: return "Anti-neutron";
        default:    return "PDG " + std::to_string(pdg);
    }
}

//********************************************************************
void pdEventDisplay::SetupCanvas(TCanvas* canvas){
//********************************************************************

    // Set up canvas style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    // Set margins
    canvas->SetLeftMargin(0.08);
    canvas->SetRightMargin(0.02);
    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.12);
}

//********************************************************************
void pdEventDisplay::DrawLegend(TLegend* legend){
//********************************************************************

    legend->SetFillColor(0);
    legend->SetBorderSize(1);
    legend->SetTextSize(0.03);
    legend->Draw();
}

//********************************************************************
bool pdEventDisplay::EventContainsRequiredParticles(AnaEventB& event) const {
//********************************************************************

    // If no required particles are specified, accept all events
    if (_RequiredParticlePDGs.empty()) {
        return true;
    }

    // Access all true particles in the event
    AnaTrueParticleB** allTrueParticles = event.TrueParticles;
    Int_t nTrueParticles = event.nTrueParticles;

    // Check if any of the required particle types are present
    for (int requiredPDG : _RequiredParticlePDGs) {
        bool foundParticle = false;

        for (int i = 0; i < nTrueParticles; i++) {
            AnaTrueParticleB* truePart = allTrueParticles[i];
            if (!truePart) continue;

            if (truePart->PDG == requiredPDG) {
                foundParticle = true;
                break;
            }
        }

        // If this required particle type is not found, reject the event
        if (!foundParticle) {
            return false;
        }
    }

    // All required particle types were found
    return true;
}

//********************************************************************
std::string pdEventDisplay::ConvertProcessToString(AnaTrueParticleB::ProcessEnum process) {
//********************************************************************
    switch (process) {
        case AnaTrueParticleB::Unknown: return "Unknown";
        case AnaTrueParticleB::primary: return "Primary";
        case AnaTrueParticleB::Decay: return "Decay";
        case AnaTrueParticleB::kaonplusInelastic: return "K+ Inelastic";
        case AnaTrueParticleB::kaonminusInelastic: return "K- Inelastic";
        case AnaTrueParticleB::neutronInelastic: return "n Inelastic";
        case AnaTrueParticleB::hadElastic: return "Hadronic Elastic";
        case AnaTrueParticleB::nCapture: return "n Capture";
        case AnaTrueParticleB::protonInelastic: return "p Inelastic";
        case AnaTrueParticleB::piplusInelastic: return "π+ Inelastic";
        case AnaTrueParticleB::piminusInelastic: return "π- Inelastic";
        case AnaTrueParticleB::hBertiniCaptureAtRest: return "Capture at Rest";
        case AnaTrueParticleB::CoulombScat: return "Coulomb Scatter";
        case AnaTrueParticleB::kaon0LInelastic: return "K0L Inelastic";
        case AnaTrueParticleB::electronNuclear: return "e Nuclear";
        case AnaTrueParticleB::muMinusCaptureAtRest: return "μ- Capture";
        case AnaTrueParticleB::dInelastic: return "d Inelastic";
        case AnaTrueParticleB::kaon0SInelastic: return "K0S Inelastic";
        case AnaTrueParticleB::positronNuclear: return "e+ Nuclear";
        case AnaTrueParticleB::lambdaInelastic: return "Λ Inelastic";
        case AnaTrueParticleB::tInelastic: return "t Inelastic";
        default: return "Other";
    }
}
