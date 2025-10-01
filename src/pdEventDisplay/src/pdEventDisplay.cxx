#include "pdEventDisplay.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TDirectory.h"
#include "TGraph2D.h"
#include "TH3F.h"
#include "ToyBoxNeutralKaon.hxx"
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
void pdEventDisplay::CreateEventDisplay(AnaEventB& event, int eventNumber, const ToyBoxNeutralKaon& box){
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
    DrawEventProjections(event, eventNumber, box);
    DrawEvent3D(event, eventNumber, box);
}

//********************************************************************
void pdEventDisplay::DrawEventProjections(AnaEventB& event, int eventNumber, const ToyBoxNeutralKaon& box){
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
    // AnaTrueParticleB* beamTrue = nullptr;
    // if (event.Beam) {
    //     AnaBeam* beam = static_cast<AnaBeam*>(event.Beam);
    //     if (beam && beam->BeamParticle) {
    //         beamTrue = beam->BeamParticle->GetTrueParticle();
    //     }
    // }

    // Create vertices using the common utility function
    // Note: maxVertexRadius parameter is not used in CreateReconstructedVertices
    // The actual vertex creation is controlled by maxDaughterDistance parameter
    const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
    std::vector<AnaVertexPD*> reconVertices = pdAnaUtils::CreateReconstructedVertices(event, _VertexRadius, maxDaughterDistance);


    // Convert AnaVertexPD* to DisplayVertex for the event display
    std::vector<DisplayVertex> vertices;
    for (const auto& reconVertex : reconVertices) {
        if (!reconVertex || reconVertex->NParticles < 2) continue; // Skip if less than 2 particles

        DisplayVertex displayVertex;

        // Set vertex position
        displayVertex.position[0] = reconVertex->Position[0];
        displayVertex.position[1] = reconVertex->Position[1];
        displayVertex.position[2] = reconVertex->Position[2];

        // Set vertex properties - these are daughter particles, not parent
        displayVertex.parent = nullptr; // No parent particle for reconstructed vertices
        displayVertex.nParticles = reconVertex->NParticles; // Number of daughter particles

        // Convert particles vector - these are the daughter particles
        displayVertex.particles.clear();
        for (int i = 0; i < reconVertex->NParticles; i++) {
            displayVertex.particles.push_back(static_cast<AnaParticlePD*>(reconVertex->Particles[i]));
        }

        displayVertex.isTrueVertex = true; // For display purposes
        displayVertex.vertexType = 1; // True vertex

        vertices.push_back(displayVertex);
    }


    // Debug: Check how many vertices have valid positions
    // int validVertices = 0;
    // for (const auto& vertex : vertices) {
    //     if (vertex.position[0] > -900 && vertex.position[1] > -900 && vertex.position[2] > -900) {
    //         validVertices++;
    //     }
    // }

    // Check if event contains signal vertices
    // Store signal flags in an array to avoid overwriting
    std::vector<bool> vertexSignalFlags(reconVertices.size(), false);


    // Create canvas with two subplots (XY and XZ projections) using event-specific name
    std::string canvasName = "c" + std::to_string(eventNumber);
    std::string canvasTitle = "Event Display - Event " + std::to_string(eventNumber) + " (Vertices: " + std::to_string(vertices.size()) + ")";

    // Add neutral particle information to the title
    if (box.neutralParticleCandidates.size() > 0) {
        canvasTitle += " [Neutral Particles: " + std::to_string(box.neutralParticleCandidates.size());

        // Add true PDG information for neutral particles
        std::string truePDGs = "";
        for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
            if (box.neutralParticleCandidates[i]->TrueObject != nullptr) {
                AnaTrueParticleB* trueNeutralParticle = static_cast<AnaTrueParticleB*>(box.neutralParticleCandidates[i]->TrueObject);
                if (!truePDGs.empty()) truePDGs += ",";
                truePDGs += std::to_string(trueNeutralParticle->PDG);
            }
        }
        if (!truePDGs.empty()) {
            canvasTitle += " (True PDGs: " + truePDGs + ")";
        }
        canvasTitle += "]";
    } else {
        canvasTitle += " [No Neutral Particles]";
    }

    TCanvas* c1 = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 1200, 600);
    c1->Divide(2, 1);

    SetupCanvas(c1);


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
    vertexGraphXY->SetMarkerStyle(21); // Square marker
    vertexGraphXY->SetMarkerSize(1.5);
    vertexGraphXY->SetMarkerColor(kBlue);

    TGraph* vertexGraphXZ = new TGraph();
    vertexGraphXZ->SetMarkerStyle(21); // Square marker
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

    // Draw circles at start position of each particle with particle color
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue;

        // Get particle color
        int particleColor = GetParticleColor(truePart->PDG);

        // Get start position using DefinePosition
        TVector3 startPos = pdAnaUtils::DefinePosition(particle);

        if (startPos.X() > -900 && startPos.Y() > -900 && startPos.Z() > -900) {
            // Draw circle in XY projection
            TEllipse* startCircleXY = new TEllipse(startPos.X(), startPos.Y(), 5.0, 5.0); // 5cm radius
            startCircleXY->SetFillStyle(0); // Hollow circle
            startCircleXY->SetLineColor(particleColor);
            startCircleXY->SetLineWidth(2);
            startCircleXY->Draw();
            vertexCirclesXY.push_back(startCircleXY); // Reuse vertex circles vector for cleanup

            // Draw circle in XZ projection
            TEllipse* startCircleXZ = new TEllipse(startPos.X(), startPos.Z(), 5.0, 5.0); // 5cm radius
            startCircleXZ->SetFillStyle(0); // Hollow circle
            startCircleXZ->SetLineColor(particleColor);
            startCircleXZ->SetLineWidth(2);
            startCircleXZ->Draw();
            vertexCirclesXZ.push_back(startCircleXZ); // Reuse vertex circles vector for cleanup

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
    std::string histTitleXY = "Event " + std::to_string(eventNumber) + " - XY Projection";

    // Add neutral particle info to histogram title
    if (box.neutralParticleCandidates.size() > 0) {
        histTitleXY += " [Neutral: " + std::to_string(box.neutralParticleCandidates.size()) + "]";
    } else {
        histTitleXY += " [No Neutral]";
    }

    histTitleXY += "; X [cm]; Y [cm]";
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

    // Draw extrapolated lines for vertices in XY projection using stored fitted parameters
    for (const auto& vertex : vertices) {
        if (vertex.position[0] < -900 || vertex.position[1] < -900) continue; // Skip invalid vertices

        // Find the matching reconstructed vertex and use its stored fitted line parameters
        for (const auto& reconVertex : reconVertices) {
            if (reconVertex->Position[0] > -900 &&
                fabs(reconVertex->Position[0] - vertex.position[0]) < 1.0 &&
                fabs(reconVertex->Position[1] - vertex.position[1]) < 1.0) {

                // Use the stored fitted line parameters instead of refitting
                for (size_t i = 0; i < reconVertex->FittedLineParams.size() && i < reconVertex->NParticles; i++) {
                    const std::vector<double>& lineParams = reconVertex->FittedLineParams[i];

                    if (lineParams.size() >= 6 && lineParams[0] != -999.0) { // Valid line parameters
                        // Get the daughter particle to determine its color
                        AnaParticlePD* daughter = static_cast<AnaParticlePD*>(reconVertex->Particles[i]);
                        int particleColor = kGreen; // Default color

                        if (daughter) {
                            // Get the true particle associated with this reconstructed particle
                            AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
                            if (trueDaughter) {
                                particleColor = GetParticleColor(trueDaughter->PDG);
                            }
                        }

                        // Draw extrapolated line in XY projection
                        // Line: P(t) = (x0, y0) + t * (dx, dy)
                        double x0 = lineParams[0];
                        double y0 = lineParams[1];
                        double dx = lineParams[3];
                        double dy = lineParams[4];

                        // Create line points extending from -200cm to +200cm in the line direction
                        double t1 = -200.0;
                        double t2 = 200.0;
                        double x1 = x0 + t1 * dx;
                        double y1 = y0 + t1 * dy;
                        double x2 = x0 + t2 * dx;
                        double y2 = y0 + t2 * dy;

                        TLine* extrapLine = new TLine(x1, y1, x2, y2);
                        extrapLine->SetLineColor(particleColor);
                        extrapLine->SetLineWidth(2);
                        extrapLine->SetLineStyle(2); // Dashed line
                        extrapLine->Draw("SAME");
                    }
                }
                break; // Found the matching vertex, no need to continue
            }
        }
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

    // Draw vertex circles using the correct radius in XY projection
    // int drawnVertices = 0;
    for (const auto& vertex : vertices) {
        // Only draw vertices with valid positions (successful line fitting)
        if (vertex.position[0] < -900 || vertex.position[1] < -900 || vertex.position[2] < -900) {
            continue; // Skip vertices with failed line fitting
        }

        // Use the calculated vertex position (where extrapolated lines intersect)
        double circleX = vertex.position[0];
        double circleY = vertex.position[1];
        int circleColor = kBlue; // Default color

        if (!vertex.particles.empty() && vertex.particles[0]) {
            // Get the color of the first daughter particle
            AnaTrueParticleB* trueDaughter = vertex.particles[0]->GetTrueParticle();
            if (trueDaughter) {
                circleColor = GetParticleColor(trueDaughter->PDG);
            }
        }


        // Draw vertex circle from particle start position with particle color
        TEllipse* vertexCircleXY = new TEllipse(circleX, circleY, maxDaughterDistance, maxDaughterDistance);
        vertexCircleXY->SetFillStyle(0); // Hollow circle
        vertexCircleXY->SetLineColor(circleColor);
        vertexCircleXY->SetLineWidth(2);
        vertexCircleXY->Draw();
        vertexCirclesXY.push_back(vertexCircleXY);
        // drawnVertices++;

    }

    // Draw particle ID labels in XY projection
    for (auto* label : particleLabelsXY) {
        if (label) label->Draw("SAME");
    }




    // Draw XZ projection (right subplot)
    c1->cd(2);
    std::string histTitleXZ = "Event " + std::to_string(eventNumber) + " - XZ Projection";

    // Add neutral particle info to histogram title
    if (box.neutralParticleCandidates.size() > 0) {
        histTitleXZ += " [Neutral: " + std::to_string(box.neutralParticleCandidates.size()) + "]";
                } else {
        histTitleXZ += " [No Neutral]";
    }

    histTitleXZ += "; X [cm]; Z [cm]";
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

    // Draw extrapolated lines for vertices in XZ projection using stored fitted parameters
    for (const auto& vertex : vertices) {
        if (vertex.position[0] < -900 || vertex.position[2] < -900) continue; // Skip invalid vertices

        // Find the matching reconstructed vertex and use its stored fitted line parameters
        for (const auto& reconVertex : reconVertices) {
            if (reconVertex->Position[0] > -900 &&
                fabs(reconVertex->Position[0] - vertex.position[0]) < 1.0 &&
                fabs(reconVertex->Position[2] - vertex.position[2]) < 1.0) {

                // Use the stored fitted line parameters instead of refitting
                for (size_t i = 0; i < reconVertex->FittedLineParams.size() && i < reconVertex->NParticles; i++) {
                    const std::vector<double>& lineParams = reconVertex->FittedLineParams[i];

                    if (lineParams.size() >= 6 && lineParams[0] != -999.0) { // Valid line parameters
                        // Get the daughter particle to determine its color
                        AnaParticlePD* daughter = static_cast<AnaParticlePD*>(reconVertex->Particles[i]);
                        int particleColor = kGreen; // Default color

                        if (daughter) {
                            // Get the true particle associated with this reconstructed particle
                            AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
                            if (trueDaughter) {
                                particleColor = GetParticleColor(trueDaughter->PDG);
                            }
                        }

                        // Draw extrapolated line in XZ projection
                        // Line: P(t) = (x0, z0) + t * (dx, dz)
                        double x0 = lineParams[0];
                        double z0 = lineParams[2];
                        double dx = lineParams[3];
                        double dz = lineParams[5];

                        // Create line points extending from -200cm to +200cm in the line direction
                        double t1 = -200.0;
                        double t2 = 200.0;
                        double x1 = x0 + t1 * dx;
                        double z1 = z0 + t1 * dz;
                        double x2 = x0 + t2 * dx;
                        double z2 = z0 + t2 * dz;

                        TLine* extrapLine = new TLine(x1, z1, x2, z2);
                        extrapLine->SetLineColor(particleColor);
                        extrapLine->SetLineWidth(2);
                        extrapLine->SetLineStyle(2); // Dashed line
                        extrapLine->Draw("SAME");
                    }
                }
                break; // Found the matching vertex, no need to continue
            }
        }
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

    // Draw vertex circles using the correct radius in XZ projection
    // int drawnVerticesXZ = 0;
    for (const auto& vertex : vertices) {
        // Only draw vertices with valid positions (successful line fitting)
        if (vertex.position[0] < -900 || vertex.position[1] < -900 || vertex.position[2] < -900) {
            continue; // Skip vertices with failed line fitting
        }

        // Use the calculated vertex position (where extrapolated lines intersect)
        double circleX = vertex.position[0];
        double circleZ = vertex.position[2];
        int circleColor = kBlue; // Default color

        if (!vertex.particles.empty() && vertex.particles[0]) {
            // Get the color of the first daughter particle
            AnaTrueParticleB* trueDaughter = vertex.particles[0]->GetTrueParticle();
            if (trueDaughter) {
                circleColor = GetParticleColor(trueDaughter->PDG);
            }
        }


        // Draw vertex circle from particle start position with particle color
        TEllipse* vertexCircleXZ = new TEllipse(circleX, circleZ, maxDaughterDistance, maxDaughterDistance);
        vertexCircleXZ->SetFillStyle(0); // Hollow circle
        vertexCircleXZ->SetLineColor(circleColor);
        vertexCircleXZ->SetLineWidth(2);
        vertexCircleXZ->Draw();
        vertexCirclesXZ.push_back(vertexCircleXZ);
        // drawnVerticesXZ++;

    }

    // Draw particle ID labels in XZ projection
    for (auto* label : particleLabelsXZ) {
        if (label) label->Draw("SAME");
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
    if (parentEndGraphXY->GetN() > 0) {
        legend->AddEntry(parentEndGraphXY, "Parent End Pos", "P");
    }
    if (parentLastHitGraphXY->GetN() > 0) {
        legend->AddEntry(parentLastHitGraphXY, "Parent Last Hit", "P");
    }
    if (daughterFirstHitGraphXY->GetN() > 0) {
        legend->AddEntry(daughterFirstHitGraphXY, "Daughter First Hit", "P");
    }

    DrawLegend(legend);








    // Draw neutral particles from the ToyBox
    const std::vector<AnaNeutralParticlePD*>& neutralParticles = box.neutralParticleCandidates;
    std::vector<TLine*> neutralParticleLinesXY;
    std::vector<TLine*> neutralParticleLinesXZ;
    std::vector<TLine*> neutralParentLinesXY;
    std::vector<TLine*> neutralParentLinesXZ;
    std::vector<TMarker*> neutralStartMarkersXY;
    std::vector<TMarker*> neutralStartMarkersXZ;
    std::vector<TMarker*> neutralEndMarkersXY;
    std::vector<TMarker*> neutralEndMarkersXZ;

    // First, draw neutral particles in XY projection
    c1->cd(1);
    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        // Get the true particle associated with the neutral particle itself
        AnaTrueParticleB* trueParticle = nullptr;

        // Check if the neutral particle has an associated true object
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        // Use the particle color based on the true PDG code
        int particleColor = GetParticleColor(trueParticle->PDG);

        // Get true start and end positions
        TVector3 startPos(trueParticle->Position[0], trueParticle->Position[1], trueParticle->Position[2]);
        TVector3 endPos(trueParticle->PositionEnd[0], trueParticle->PositionEnd[1], trueParticle->PositionEnd[2]);

        // Check if positions are valid
        if (startPos.X() < -900 || endPos.X() < -900) continue;

        // Draw line in XY projection
        TLine* lineXY = new TLine(startPos.X(), startPos.Y(), endPos.X(), endPos.Y());
        lineXY->SetLineColor(particleColor);
        lineXY->SetLineWidth(2);
        lineXY->SetLineStyle(1); // Solid line
        lineXY->Draw();
        neutralParticleLinesXY.push_back(lineXY);

        // Draw start marker (circle) in XY projection
        TMarker* startMarkerXY = new TMarker(startPos.X(), startPos.Y(), 20); // Circle marker
        startMarkerXY->SetMarkerColor(particleColor);
        startMarkerXY->SetMarkerSize(1.2);
        startMarkerXY->Draw();
        neutralStartMarkersXY.push_back(startMarkerXY);

        // Draw end marker (square) in XY projection
        TMarker* endMarkerXY = new TMarker(endPos.X(), endPos.Y(), 21); // Square marker
        endMarkerXY->SetMarkerColor(particleColor);
        endMarkerXY->SetMarkerSize(1.2);
        endMarkerXY->Draw();
        neutralEndMarkersXY.push_back(endMarkerXY);

        // Draw parent particle fitted line if it exists
        AnaParticlePD* parentParticle = neutralParticle->Parent;
        if (parentParticle) {
            AnaTrueParticleB* parentTrueParticle = parentParticle->GetTrueParticle();
            if (parentTrueParticle) {
                // Use the same color as the neutral particle for consistency
                int parentColor = particleColor;

                // Get parent particle end position (for parent particles, use end position)
                TVector3 parentStartPos = pdAnaUtils::DefinePosition(parentParticle, false);

                // Get parent particle direction (from momentum)
                TVector3 parentDirection(parentTrueParticle->Direction[0],
                                        parentTrueParticle->Direction[1],
                                        parentTrueParticle->Direction[2]);

                // Create parent line extending from start position
                double lineLength = 100.0; // cm
                TVector3 parentEndPos = parentStartPos + lineLength * parentDirection;

                // Draw parent line in XY projection
                TLine* parentLineXY = new TLine(parentStartPos.X(), parentStartPos.Y(),
                                              parentEndPos.X(), parentEndPos.Y());
                parentLineXY->SetLineColor(parentColor);
                parentLineXY->SetLineWidth(2);
                parentLineXY->SetLineStyle(3); // Dotted line to distinguish from neutral particle
                parentLineXY->Draw();
                neutralParentLinesXY.push_back(parentLineXY);
            }
        }
    }

    // Then, draw neutral particles in XZ projection
    c1->cd(2);
    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        // Get the true particle associated with the neutral particle itself
        AnaTrueParticleB* trueParticle = nullptr;

        // Check if the neutral particle has an associated true object
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        // Use the particle color based on the true PDG code
        int particleColor = GetParticleColor(trueParticle->PDG);

        // Get true start and end positions
        TVector3 startPos(trueParticle->Position[0], trueParticle->Position[1], trueParticle->Position[2]);
        TVector3 endPos(trueParticle->PositionEnd[0], trueParticle->PositionEnd[1], trueParticle->PositionEnd[2]);

        // Check if positions are valid
        if (startPos.X() < -900 || endPos.X() < -900) continue;

        // Draw line in XZ projection
        TLine* lineXZ = new TLine(startPos.X(), startPos.Z(), endPos.X(), endPos.Z());
        lineXZ->SetLineColor(particleColor);
        lineXZ->SetLineWidth(2);
        lineXZ->SetLineStyle(1); // Solid line
        lineXZ->Draw();
        neutralParticleLinesXZ.push_back(lineXZ);

        // Draw start marker (circle) in XZ projection
        TMarker* startMarkerXZ = new TMarker(startPos.X(), startPos.Z(), 20); // Circle marker
        startMarkerXZ->SetMarkerColor(particleColor);
        startMarkerXZ->SetMarkerSize(1.2);
        startMarkerXZ->Draw();
        neutralStartMarkersXZ.push_back(startMarkerXZ);

        // Draw end marker (square) in XZ projection
        TMarker* endMarkerXZ = new TMarker(endPos.X(), endPos.Z(), 21); // Square marker
        endMarkerXZ->SetMarkerColor(particleColor);
        endMarkerXZ->SetMarkerSize(1.2);
        endMarkerXZ->Draw();
        neutralEndMarkersXZ.push_back(endMarkerXZ);

        // Draw parent particle fitted line if it exists
        AnaParticlePD* parentParticle = neutralParticle->Parent;
        if (parentParticle) {
            AnaTrueParticleB* parentTrueParticle = parentParticle->GetTrueParticle();
            if (parentTrueParticle) {
                // Use the same color as the neutral particle for consistency
                int parentColor = particleColor;

                // Get parent particle end position (for parent particles, use end position)
                TVector3 parentStartPos = pdAnaUtils::DefinePosition(parentParticle, false);

                // Get parent particle direction (from momentum)
                TVector3 parentDirection(parentTrueParticle->Direction[0],
                                        parentTrueParticle->Direction[1],
                                        parentTrueParticle->Direction[2]);

                // Create parent line extending from start position
                double lineLength = 100.0; // cm
                TVector3 parentEndPos = parentStartPos + lineLength * parentDirection;

                // Draw parent line in XZ projection
                TLine* parentLineXZ = new TLine(parentStartPos.X(), parentStartPos.Z(),
                                              parentEndPos.X(), parentEndPos.Z());
                parentLineXZ->SetLineColor(parentColor);
                parentLineXZ->SetLineWidth(2);
                parentLineXZ->SetLineStyle(3); // Dotted line to distinguish from neutral particle
                parentLineXZ->Draw();
                neutralParentLinesXZ.push_back(parentLineXZ);
            }
        }
    }

    // Add individual legend entries for neutral particle parent lines
    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        // Get the true particle associated with the neutral particle itself
        AnaTrueParticleB* trueParticle = nullptr;
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        // Use the particle color based on the true PDG code
        int particleColor = GetParticleColor(trueParticle->PDG);

        // Create legend entry for parent line
        TLine* dummyParentLine = new TLine();
        dummyParentLine->SetLineStyle(3); // Dotted line
        dummyParentLine->SetLineColor(particleColor);
        dummyParentLine->SetLineWidth(2);

        // Create legend text with particle information
        std::string legendText = "Parent #" + std::to_string(i+1) + " (PDG: " + std::to_string(trueParticle->PDG) + ")";
        legend->AddEntry(dummyParentLine, legendText.c_str(), "L");
    }

    // Redraw legend with all parent entries
    if (!neutralParentLinesXY.empty()) {
        legend->Draw();
    }

    // Add individual legend entries for each neutral particle
    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        // Get the true particle associated with the neutral particle itself
        AnaTrueParticleB* trueParticle = nullptr;
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        // Use the particle color based on the true PDG code
        int particleColor = GetParticleColor(trueParticle->PDG);

        // Create legend entry with particle information
        TLine* dummyNeutralLine = new TLine();
        dummyNeutralLine->SetLineStyle(1); // Solid line
        dummyNeutralLine->SetLineColor(particleColor);
        dummyNeutralLine->SetLineWidth(2);

        // Create legend text with particle PDG and index
        std::string legendText = "Neutral Particle #" + std::to_string(i+1) + " (PDG: " + std::to_string(trueParticle->PDG) + ")";
        legend->AddEntry(dummyNeutralLine, legendText.c_str(), "L");
    }

    // Redraw legend with all neutral particle entries
    if (!neutralParticleLinesXY.empty()) {
        legend->Draw();
    }

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


    // Clean up neutral particle objects
    for (auto* line : neutralParticleLinesXY) {
        delete line;
    }
    for (auto* line : neutralParticleLinesXZ) {
        delete line;
    }
    for (auto* line : neutralParentLinesXY) {
        delete line;
    }
    for (auto* line : neutralParentLinesXZ) {
        delete line;
    }
    for (auto* marker : neutralStartMarkersXY) {
        delete marker;
    }
    for (auto* marker : neutralStartMarkersXZ) {
        delete marker;
    }
    for (auto* marker : neutralEndMarkersXY) {
        delete marker;
    }
    for (auto* marker : neutralEndMarkersXZ) {
        delete marker;
    }

    // Clean up particle labels
    for (auto* label : particleLabelsXY) {
        if (label) delete label;
    }
    for (auto* label : particleLabelsXZ) {
        if (label) delete label;
    }


    // Note: We don't delete the canvas here as it might be displayed on screen
    // The canvas will be cleaned up by ROOT's garbage collection
}

//********************************************************************
void pdEventDisplay::DrawEvent3D(AnaEventB& event, int eventNumber, const ToyBoxNeutralKaon& box){
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

    // Create vertices using the common utility function
    const double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
    std::vector<AnaVertexPD*> reconVertices = pdAnaUtils::CreateReconstructedVertices(event, _VertexRadius, maxDaughterDistance);

    // Convert reconstructed vertices to display vertices
    std::vector<DisplayVertex> vertices;
    for (const auto& reconVertex : reconVertices) {
        if (!reconVertex) continue;

        DisplayVertex displayVertex;
        displayVertex.position[0] = reconVertex->Position[0];
        displayVertex.position[1] = reconVertex->Position[1];
        displayVertex.position[2] = reconVertex->Position[2];
        displayVertex.parent = static_cast<AnaParticlePD*>(reconVertex->Particles[0]);
        displayVertex.nParticles = reconVertex->NParticles;
        displayVertex.isTrueVertex = false;
        displayVertex.vertexType = 0; // reconstructed

        // Add all particles associated with this vertex
        for (int i = 0; i < reconVertex->NParticles; i++) {
            AnaParticlePD* vertexParticle = static_cast<AnaParticlePD*>(reconVertex->Particles[i]);
            if (vertexParticle) {
                displayVertex.particles.push_back(vertexParticle);
            }
        }

        vertices.push_back(displayVertex);
    }

    // Create canvas with 3D view
    std::string canvasName = "c3d_" + std::to_string(eventNumber);
    std::string canvasTitle = "3D Event Display - Event " + std::to_string(eventNumber) + " (Vertices: " + std::to_string(vertices.size()) + ") - Left: rotate, Right: zoom, Middle: pan";

    // Add neutral particle information to the title
    if (box.neutralParticleCandidates.size() > 0) {
        canvasTitle += " [Neutral Particles: " + std::to_string(box.neutralParticleCandidates.size());

        // Add true PDG information for neutral particles
        std::string truePDGs = "";
        for (size_t i = 0; i < box.neutralParticleCandidates.size(); i++) {
            if (box.neutralParticleCandidates[i]->TrueObject != nullptr) {
                AnaTrueParticleB* trueNeutralParticle = static_cast<AnaTrueParticleB*>(box.neutralParticleCandidates[i]->TrueObject);
                if (!truePDGs.empty()) truePDGs += ",";
                truePDGs += std::to_string(trueNeutralParticle->PDG);
            }
        }
        if (!truePDGs.empty()) {
            canvasTitle += " (True PDGs: " + truePDGs + ")";
        }
        canvasTitle += "]";
    } else {
        canvasTitle += " [No Neutral Particles]";
    }

    TCanvas* c3d = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 600);
    SetupCanvas(c3d);

    // Create 3D histogram for coordinate system
    TH3F* h3d = new TH3F("h3d", "3D Event Display", 50, -360, 360, 50, 0, 700, 50, 0, 700);
    h3d->SetStats(0);
    h3d->GetXaxis()->SetTitle("X [cm]");
    h3d->GetYaxis()->SetTitle("Y [cm]");
    h3d->GetZaxis()->SetTitle("Z [cm]");
    h3d->Draw();

    // Set up 3D view with proper controls
    c3d->SetPhi(30);  // Azimuthal angle (0-360 degrees)
    c3d->SetTheta(30); // Polar angle (0-180 degrees)

    // Enable 3D interaction
    c3d->SetEditable(kTRUE);
    c3d->ToggleEventStatus();

    // Set the view to 3D mode
    TView3D* view3d = (TView3D*)TView::CreateView(1);
    view3d->SetRange(-360, 0, 0, 360, 700, 700);
    Int_t irep = 0;
    view3d->SetView(0.5, 0.5, 0.5, irep);
    c3d->SetView(view3d);

    // Configure canvas for better 3D interaction
    c3d->SetCanvasSize(800, 600);
    c3d->SetRightMargin(0.05);
    c3d->SetLeftMargin(0.05);
    c3d->SetTopMargin(0.05);
    c3d->SetBottomMargin(0.05);

    // Configure 3D histogram ranges - these can be manually adjusted
    h3d->GetXaxis()->SetRangeUser(-360, 360);
    h3d->GetYaxis()->SetRangeUser(0, 700);
    h3d->GetZaxis()->SetRangeUser(0, 700);

    // Set up for manual zooming by setting specific ranges
    // You can modify these values to zoom in on specific regions
    // Example: to zoom in on X-axis from -100 to 100:
    // h3d->GetXaxis()->SetRangeUser(-100, 100);

    // Add instructions to canvas title for manual range setting
    canvasTitle += " - Use Set3DAxisRanges(h3d, xmin, xmax, ymin, ymax, zmin, zmax) to zoom";

    // Example zoom presets (uncomment to use):
    // Set3DAxisRanges(h3d, -100, 100, 200, 400, 200, 400);  // Zoom to center region
    // Set3DAxisRanges(h3d, -200, 200, 0, 200, 0, 200);      // Zoom to front corner
    // Set3DAxisRanges(h3d, -50, 50, 300, 500, 300, 500);    // Zoom to specific vertex region

    // Test zoom (uncomment to test):
    // Set3DAxisRanges(h3d, -50, 50, 300, 500, 300, 500);    // Test zoom to center region

    // Storage for 3D objects
    std::vector<TGraph2D*> eventParticleGraphs3D;
    std::vector<TGraph2D*> parentEndGraph3D;
    std::vector<TGraph2D*> parentLastHitGraph3D;
    std::vector<TGraph2D*> daughterFirstHitGraph3D;
    std::vector<TGraph2D*> daughterStartGraph3D;
    std::vector<TPolyLine3D*> vertexLines3D;
    std::vector<TGraph2D*> vertexMarkers3D;
    std::vector<TPolyLine3D*> neutralParticleLines3D;
    std::vector<TPolyLine3D*> neutralParentLines3D;
    std::vector<TGraph2D*> neutralStartMarkers3D;
    std::vector<TGraph2D*> neutralEndMarkers3D;

    // Create graphs for different particle types
    std::map<int, TGraph2D*> eventParticleGraphsMap3D;

    // Fill 3D graphs with hit data
    for (const auto& particle : allParticles) {
        AnaTrueParticleB* truePart = particle->GetTrueParticle();
        if (!truePart) continue;

        int particleType = truePart->PDG;
        if (eventParticleGraphsMap3D.find(particleType) == eventParticleGraphsMap3D.end()) {
            int color = GetParticleColor(particleType);
            TGraph2D* graph3D = new TGraph2D();
            graph3D->SetMarkerStyle(20);
            graph3D->SetMarkerSize(0.5);
            graph3D->SetMarkerColor(color);
            eventParticleGraphsMap3D[particleType] = graph3D;
        }

        // Add hits to 3D graph
        for (int plane = 0; plane < 3; plane++) {
            for (UInt_t h = 0; h < particle->Hits[plane].size(); h++) {
                AnaHitPD hit = particle->Hits[plane][h];
                if (hit.Position.X() > -900 && hit.Position.Y() > -900 && hit.Position.Z() > -900) {
                    int nPoints = eventParticleGraphsMap3D[particleType]->GetN();
                    eventParticleGraphsMap3D[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Y(), hit.Position.Z());
                }
            }
        }
    }

    // Draw all particle graphs
    for (auto& pair : eventParticleGraphsMap3D) {
        if (pair.second->GetN() > 0) {
            pair.second->Draw("P SAME");
        }
    }

    // Draw extrapolated tracks for vertices
    for (size_t v = 0; v < vertices.size(); v++) {
        const auto& vertex = vertices[v];
        if (vertex.position[0] < -900 || vertex.position[1] < -900 || vertex.position[2] < -900) {
            continue; // Skip vertices with failed line fitting
        }

        // Get the corresponding reconstructed vertex to access FittedLineParams
        if (v < reconVertices.size() && reconVertices[v]) {
            AnaVertexPD* reconVertex = reconVertices[v];

            // Draw extrapolated tracks for each daughter particle
            for (size_t i = 0; i < vertex.particles.size() && i < reconVertex->FittedLineParams.size(); i++) {
                const auto& daughter = vertex.particles[i];
                if (!daughter) continue;

                // Get the fitted line parameters from the vertex
                const auto& lineParams = reconVertex->FittedLineParams[i];
                if (lineParams.size() >= 6) {
                    double x0 = lineParams[0];
                    double y0 = lineParams[1];
                    double z0 = lineParams[2];
                    double dx = lineParams[3];
                    double dy = lineParams[4];
                    double dz = lineParams[5];

                    // Get the true particle for color
                    AnaTrueParticleB* truePart = daughter->GetTrueParticle();
                    int particleColor = kBlue; // Default color
                    if (truePart) {
                        particleColor = GetParticleColor(truePart->PDG);
                    }

                    // Create extrapolated track line
                    TPolyLine3D* extrapolatedTrack = new TPolyLine3D(2);

                    // Calculate start and end points of the extrapolated line across the full detector range
                    // Detector range: x[-360,360], y[0,700], z[0,700]
                    double xMin = -360, xMax = 360;
                    double yMin = 0, yMax = 700;
                    double zMin = 0, zMax = 700;

                    // Calculate parameter t for each boundary
                    double tXMin = (xMin - x0) / dx;
                    double tXMax = (xMax - x0) / dx;
                    double tYMin = (yMin - y0) / dy;
                    double tYMax = (yMax - y0) / dy;
                    double tZMin = (zMin - z0) / dz;
                    double tZMax = (zMax - z0) / dz;

                    // Find the valid range of t that keeps the line within detector bounds
                    double tMin = std::max({tXMin, tYMin, tZMin});
                    double tMax = std::min({tXMax, tYMax, tZMax});

                    // Ensure we have a valid range
                    if (tMin < tMax) {
                        double startX = x0 + dx * tMin;
                        double startY = y0 + dy * tMin;
                        double startZ = z0 + dz * tMin;
                        double endX = x0 + dx * tMax;
                        double endY = y0 + dy * tMax;
                        double endZ = z0 + dz * tMax;

                        extrapolatedTrack->SetPoint(0, startX, startY, startZ);
                        extrapolatedTrack->SetPoint(1, endX, endY, endZ);
                        extrapolatedTrack->SetLineColor(particleColor);
                        extrapolatedTrack->SetLineWidth(2);
                        extrapolatedTrack->SetLineStyle(2); // Dashed line
                        extrapolatedTrack->Draw();
                        vertexLines3D.push_back(extrapolatedTrack);
                    }
                }
            }
        }
    }

    // Draw vertices
    for (const auto& vertex : vertices) {
        if (vertex.position[0] < -900 || vertex.position[1] < -900 || vertex.position[2] < -900) {
            continue; // Skip vertices with failed line fitting
        }

        // Draw vertex position marker
        TGraph2D* vertexMarker = new TGraph2D();
        vertexMarker->SetPoint(0, vertex.position[0], vertex.position[1], vertex.position[2]);
        vertexMarker->SetMarkerStyle(kFullSquare);
        vertexMarker->SetMarkerColor(kBlue);
        vertexMarker->SetMarkerSize(1.5);
        vertexMarker->Draw("P SAME");
        vertexMarkers3D.push_back(vertexMarker);

        // Draw lines from vertex to daughter particles
        for (const auto& daughter : vertex.particles) {
            if (!daughter) continue;

            // Get daughter start position
            double daughterX = daughter->PositionStart[0];
            double daughterY = daughter->PositionStart[1];
            double daughterZ = daughter->PositionStart[2];

            if (daughterX > -900 && daughterY > -900 && daughterZ > -900) {
                TPolyLine3D* vertexLine = new TPolyLine3D(2);
                vertexLine->SetPoint(0, vertex.position[0], vertex.position[1], vertex.position[2]);
                vertexLine->SetPoint(1, daughterX, daughterY, daughterZ);
                vertexLine->SetLineColor(kBlue);
                vertexLine->SetLineWidth(2);
                vertexLine->Draw();
                vertexLines3D.push_back(vertexLine);
            }
        }
    }

    // Draw neutral particles from the ToyBox
    const std::vector<AnaNeutralParticlePD*>& neutralParticles = box.neutralParticleCandidates;

    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        // Get the true particle associated with the neutral particle itself
        AnaTrueParticleB* trueParticle = nullptr;
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        // Use the particle color based on the true PDG code
        int particleColor = GetParticleColor(trueParticle->PDG);

        // Draw true neutral particle line
        TPolyLine3D* neutralLine = new TPolyLine3D(2);
        neutralLine->SetPoint(0, trueParticle->Position[0], trueParticle->Position[1], trueParticle->Position[2]);
        neutralLine->SetPoint(1, trueParticle->PositionEnd[0], trueParticle->PositionEnd[1], trueParticle->PositionEnd[2]);
        neutralLine->SetLineColor(particleColor);
        neutralLine->SetLineWidth(3);
        neutralLine->SetLineStyle(1);
        neutralLine->Draw();
        neutralParticleLines3D.push_back(neutralLine);

        // Draw start and end markers
        TGraph2D* startMarker = new TGraph2D();
        startMarker->SetPoint(0, trueParticle->Position[0], trueParticle->Position[1], trueParticle->Position[2]);
        startMarker->SetMarkerStyle(kFullCircle);
        startMarker->SetMarkerColor(particleColor);
        startMarker->SetMarkerSize(1.0);
        startMarker->Draw("P SAME");
        neutralStartMarkers3D.push_back(startMarker);

        TGraph2D* endMarker = new TGraph2D();
        endMarker->SetPoint(0, trueParticle->PositionEnd[0], trueParticle->PositionEnd[1], trueParticle->PositionEnd[2]);
        endMarker->SetMarkerStyle(kFullSquare);
        endMarker->SetMarkerColor(particleColor);
        endMarker->SetMarkerSize(1.0);
        endMarker->Draw("P SAME");
        neutralEndMarkers3D.push_back(endMarker);

        // Draw parent line if available
        if (neutralParticle->Parent) {
            AnaParticlePD* parent = neutralParticle->Parent;
            TPolyLine3D* parentLine = new TPolyLine3D(2);
            parentLine->SetPoint(0, parent->PositionStart[0], parent->PositionStart[1], parent->PositionStart[2]);
            parentLine->SetPoint(1, parent->PositionEnd[0], parent->PositionEnd[1], parent->PositionEnd[2]);
            parentLine->SetLineColor(particleColor);
            parentLine->SetLineWidth(2);
            parentLine->SetLineStyle(3); // Dotted line
            parentLine->Draw();
            neutralParentLines3D.push_back(parentLine);
        }
    }

    // Create legend
    TLegend* legend3D = new TLegend(0.65, 0.65, 0.95, 0.95);
    legend3D->SetNColumns(1);
    legend3D->SetTextSize(0.025);

    // Add legend entries for particles that have hits
    for (auto& pair : eventParticleGraphsMap3D) {
        if (pair.second->GetN() > 0) {
            std::string particleName = GetParticleName(pair.first);
            legend3D->AddEntry(pair.second, particleName.c_str(), "P");
        }
    }

    // Add neutral particle legend entries
    for (size_t i = 0; i < neutralParticles.size(); i++) {
        const auto& neutralParticle = neutralParticles[i];
        if (!neutralParticle) continue;

        AnaTrueParticleB* trueParticle = nullptr;
        if (neutralParticle->TrueObject != nullptr) {
            trueParticle = static_cast<AnaTrueParticleB*>(neutralParticle->TrueObject);
        }

        if (!trueParticle) continue;

        int particleColor = GetParticleColor(trueParticle->PDG);

        // Create legend entry for neutral particle
        TPolyLine3D* dummyNeutralLine = new TPolyLine3D(2);
        dummyNeutralLine->SetLineStyle(1);
        dummyNeutralLine->SetLineColor(particleColor);
        dummyNeutralLine->SetLineWidth(2);

        std::string legendText = "Neutral Particle #" + std::to_string(i+1) + " (PDG: " + std::to_string(trueParticle->PDG) + ")";
        legend3D->AddEntry(dummyNeutralLine, legendText.c_str(), "L");
    }

    // Add vertex legend entry
    if (vertexMarkers3D.size() > 0) {
        legend3D->AddEntry(vertexMarkers3D[0], "Vertex Position", "P");
    }

    // Add extrapolated tracks legend entry
    if (vertexLines3D.size() > 0) {
        TPolyLine3D* dummyExtrapolatedLine = new TPolyLine3D(2);
        dummyExtrapolatedLine->SetLineStyle(2); // Dashed line
        dummyExtrapolatedLine->SetLineColor(kBlack);
        dummyExtrapolatedLine->SetLineWidth(2);
        legend3D->AddEntry(dummyExtrapolatedLine, "Extrapolated Tracks", "L");
    }

    DrawLegend(legend3D);

    // Save canvas to ROOT file if enabled
    if (_SaveToRootFile) {
        SaveCanvasToRootFile(c3d, eventNumber);
    }

    // Clean up
    for (auto& pair : eventParticleGraphsMap3D) {
        if (pair.second) delete pair.second;
    }
    for (auto* line : vertexLines3D) {
        if (line) delete line;
    }
    for (auto* marker : vertexMarkers3D) {
        if (marker) delete marker;
    }
    for (auto* line : neutralParticleLines3D) {
        if (line) delete line;
    }
    for (auto* line : neutralParentLines3D) {
        if (line) delete line;
    }
    for (auto* marker : neutralStartMarkers3D) {
        if (marker) delete marker;
    }
    for (auto* marker : neutralEndMarkers3D) {
        if (marker) delete marker;
    }
    if (h3d) delete h3d;
}

// TODO: Add method implementations here

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

// Helper method to set 3D axis ranges for manual zooming
void pdEventDisplay::Set3DAxisRanges(TH3F* h3d, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    if (!h3d) return;

    // Set the axis ranges
    h3d->GetXaxis()->SetRangeUser(xmin, xmax);
    h3d->GetYaxis()->SetRangeUser(ymin, ymax);
    h3d->GetZaxis()->SetRangeUser(zmin, zmax);

    // Force the histogram to update
    h3d->SetAxisRange(xmin, xmax, "X");
    h3d->SetAxisRange(ymin, ymax, "Y");
    h3d->SetAxisRange(zmin, zmax, "Z");

    // Update the 3D view range
    TView3D* view3d = (TView3D*)gPad->GetView();
    if (view3d) {
        view3d->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
    }

    // Alternative approach: recreate the view
    if (!view3d) {
        TView3D* newView3d = (TView3D*)TView::CreateView(1);
        newView3d->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
        Int_t irep = 0;
        newView3d->SetView(0.5, 0.5, 0.5, irep);
        gPad->SetView(newView3d);
    }

    // Update the display
    h3d->Draw();
    gPad->Update();
}

