// File: LongPionsK0Length.C
// Script to calculate distances between K0 parent end position and π+/- start positions
// Uses truth information to find K0 and parent, but calculates distances using RECO information
// Only processes events with π+ and π- daughters that travel at least 10cm each
void LongPionsK0Length(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
    // Load necessary libraries
    gSystem->Load("libhighland");
    gSystem->Load("libhighlandPD");

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
        f->Close();
        return;
    }

    // Set up the tree
    AnaSpillB* spill = new AnaSpillB();
    MiniTree->SetBranchAddress("Spill", &spill);

    Long64_t nEntries = MiniTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;

    // Event counter
    int longPionEventCount = 0;

    // Create histograms
    TH1F* hPiPlusDistance = new TH1F("hPiPlusDistance", "Distance from K0 Parent to #pi^{+} Start;Distance [cm];Entries", 50, 0, 200);
    hPiPlusDistance->SetLineColor(kBlue);
    hPiPlusDistance->SetLineWidth(2);

    TH1F* hPiMinusDistance = new TH1F("hPiMinusDistance", "Distance from K0 Parent to #pi^{-} Start;Distance [cm];Entries", 50, 0, 200);
    hPiMinusDistance->SetLineColor(kRed);
    hPiMinusDistance->SetLineWidth(2);

    TH1F* hPiPlusAngle = new TH1F("hPiPlusAngle", "Angle between K0 Parent->#pi^{+} and #pi^{+} Direction;#pi^{+} Angle [degrees];Entries", 50, 0, 180);
    hPiPlusAngle->SetLineColor(kBlue);
    hPiPlusAngle->SetLineWidth(2);

    TH1F* hPiMinusAngle = new TH1F("hPiMinusAngle", "Angle between K0 Parent->#pi^{-} and #pi^{-} Direction;#pi^{-} Angle [degrees];Entries", 50, 0, 180);
    hPiMinusAngle->SetLineColor(kRed);
    hPiMinusAngle->SetLineWidth(2);

    TH1F* hCombinedAngle = new TH1F("hCombinedAngle", "Angle between K0 Direction and Pion System Direction;Angle [degrees];Entries", 50, 0, 180);
    hCombinedAngle->SetLineColor(kGreen);
    hCombinedAngle->SetLineWidth(2);

    // 2D histograms for distance vs angle relationships
    TH2F* hPiPlusDistanceVsAngle = new TH2F("hPiPlusDistanceVsAngle", "Distance vs Angle: #pi^{+};Distance [cm];Angle [degrees]", 50, 0, 200, 50, 0, 180);
    hPiPlusDistanceVsAngle->SetMarkerStyle(20);
    hPiPlusDistanceVsAngle->SetMarkerSize(0.5);
    hPiPlusDistanceVsAngle->SetMarkerColor(kBlue);

    TH2F* hPiMinusDistanceVsAngle = new TH2F("hPiMinusDistanceVsAngle", "Distance vs Angle: #pi^{-};Distance [cm];Angle [degrees]", 50, 0, 200, 50, 0, 180);
    hPiMinusDistanceVsAngle->SetMarkerStyle(20);
    hPiMinusDistanceVsAngle->SetMarkerSize(0.5);
    hPiMinusDistanceVsAngle->SetMarkerColor(kRed);

    TH2F* hPiPlusDistanceVsCombinedAngle = new TH2F("hPiPlusDistanceVsCombinedAngle", "Distance vs K0-Pion System Angle: #pi^{+};Distance [cm];Angle [degrees]", 50, 0, 200, 50, 0, 180);
    hPiPlusDistanceVsCombinedAngle->SetMarkerStyle(20);
    hPiPlusDistanceVsCombinedAngle->SetMarkerSize(0.5);
    hPiPlusDistanceVsCombinedAngle->SetMarkerColor(kBlue);

    TH2F* hPiMinusDistanceVsCombinedAngle = new TH2F("hPiMinusDistanceVsCombinedAngle", "Distance vs K0-Pion System Angle: #pi^{-};Distance [cm];Angle [degrees]", 50, 0, 200, 50, 0, 180);
    hPiMinusDistanceVsCombinedAngle->SetMarkerStyle(20);
    hPiMinusDistanceVsCombinedAngle->SetMarkerSize(0.5);
    hPiMinusDistanceVsCombinedAngle->SetMarkerColor(kRed);

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);

        // First, find all K0 particles and their daughters based on TRUTH information
        std::vector<AnaTrueParticleB*> k0Particles;
        std::vector<AnaTrueParticleB*> k0PiPlusDaughters;
        std::vector<AnaTrueParticleB*> k0PiMinusDaughters;

        // Find all K0 particles
        for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
            AnaTrueParticleB* truePart = spill->TrueParticles[j];
            if (!truePart) continue;
            if (truePart->PDG == 310) { // K0
                k0Particles.push_back(truePart);
            }
        }

        // Find daughters of each K0 and check for exactly one π+ and one π-
        for (UInt_t k0Idx = 0; k0Idx < k0Particles.size(); k0Idx++) {
            AnaTrueParticleB* k0Particle = k0Particles[k0Idx];

            // Count daughters of this K0
            int piPlusCount = 0;
            int piMinusCount = 0;
            int totalDaughterCount = 0;

            for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Check if this particle is a daughter of the K0 (using ParentID)
                if (truePart->ParentID == k0Particle->ID) {
                    totalDaughterCount++;
                    if (truePart->PDG == 211) {
                        piPlusCount++;
                        k0PiPlusDaughters.push_back(truePart);
                    } else if (truePart->PDG == -211) {
                        piMinusCount++;
                        k0PiMinusDaughters.push_back(truePart);
                    }
                }
            }

            // Only keep this K0 if it has exactly 2 daughters: 1 π+ and 1 π-
            if (totalDaughterCount == 2 && piPlusCount == 1 && piMinusCount == 1) {
                // Print K0 decay reaction
                std::cout << "Event " << i << ": K0 decay found - K0(ID:" << k0Particle->ID << ") -> π+(ID:" << k0PiPlusDaughters[0]->ID << ") + π-(ID:" << k0PiMinusDaughters[0]->ID << ")" << std::endl;

                // This K0 qualifies - keep the daughters we found
                break; // We only need one qualifying K0 per event
            } else {
                // This K0 doesn't qualify - clear the daughters we found
                k0PiPlusDaughters.clear();
                k0PiMinusDaughters.clear();
            }
        }

        // Now find reconstructed particles that match these truth daughters
        bool hasLongPiPlus = false;
        bool hasLongPiMinus = false;
        std::vector<AnaParticleB*> longPiPlusParticles;
        std::vector<AnaParticleB*> longPiMinusParticles;

        // Check reconstructed particles for long pions that are K0 daughters
        for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
            AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
            if (!bunch) continue;

            for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                AnaParticleB* particle = bunch->Particles[j];
                if (!particle) continue;

                AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                if (!particlePD) continue;

                // Get associated true particle
                AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                if (!truePart) continue;

                // Check if this reconstructed particle corresponds to a K0 daughter
                bool isK0PiPlusDaughter = false;
                bool isK0PiMinusDaughter = false;

                for (UInt_t k = 0; k < k0PiPlusDaughters.size(); k++) {
                    if (k0PiPlusDaughters[k]->ID == truePart->ID) {
                        isK0PiPlusDaughter = true;
                        break;
                    }
                }

                for (UInt_t k = 0; k < k0PiMinusDaughters.size(); k++) {
                    if (k0PiMinusDaughters[k]->ID == truePart->ID) {
                        isK0PiMinusDaughter = true;
                        break;
                    }
                }

                // If this is a K0 daughter, check track length
                if (isK0PiPlusDaughter || isK0PiMinusDaughter) {
                    // Calculate track length from start and end positions
                    float startX = particlePD->PositionStart[0];
                    float startY = particlePD->PositionStart[1];
                    float startZ = particlePD->PositionStart[2];
                    float endX = particlePD->PositionEnd[0];
                    float endY = particlePD->PositionEnd[1];
                    float endZ = particlePD->PositionEnd[2];

                    if (startX > -900 && startY > -900 && startZ > -900 &&
                        endX > -900 && endY > -900 && endZ > -900) {

                        float trackLength = sqrt((endX - startX)*(endX - startX) +
                                               (endY - startY)*(endY - startY) +
                                               (endZ - startZ)*(endZ - startZ));

                        if (trackLength > 10.0) { // 10cm threshold
                            if (isK0PiPlusDaughter) {
                                hasLongPiPlus = true;
                                longPiPlusParticles.push_back(particle);
                            } else if (isK0PiMinusDaughter) {
                                hasLongPiMinus = true;
                                longPiMinusParticles.push_back(particle);
                            }
                        }
                    }
                }
            }
        }

        // Only process events with exactly one π+ and one π- daughter of K0, both traveling >10cm
        if (k0PiPlusDaughters.size() == 1 && k0PiMinusDaughters.size() == 1 && hasLongPiPlus && hasLongPiMinus) {
            longPionEventCount++;

            // Print detailed decay information for this event
            std::cout << "=== Event " << i << " (K0 Decay Event #" << longPionEventCount << ") ===" << std::endl;

            // Find the K0 and its parent for distance calculation
            AnaTrueParticleB* k0Particle = nullptr;
            AnaTrueParticleB* k0Parent = nullptr;
            for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;
                if (truePart->PDG == 310) { // K0
                    k0Particle = truePart;
                    // Find the parent of this K0
                    for (UInt_t k = 0; k < spill->TrueParticles.size(); k++) {
                        AnaTrueParticleB* parentPart = spill->TrueParticles[k];
                        if (!parentPart) continue;
                        if (parentPart->ID == truePart->ParentID) {
                            k0Parent = parentPart;
                            break;
                        }
                    }
                    break;
                }
            }

            if (k0Particle && k0Parent) {
                // Print complete decay chain
                std::cout << "  Decay chain: Parent(ID:" << k0Parent->ID << ", PDG:" << k0Parent->PDG << ") -> K0(ID:" << k0Particle->ID << ") -> π+(ID:" << k0PiPlusDaughters[0]->ID << ") + π-(ID:" << k0PiMinusDaughters[0]->ID << ")" << std::endl;

                // Get K0 parent end position (from truth)
                float parentEndX = k0Parent->PositionEnd[0];
                float parentEndY = k0Parent->PositionEnd[1];
                float parentEndZ = k0Parent->PositionEnd[2];

                // Calculate combined direction vector from the two pions (sum of their direction vectors)
                float combinedDirX = 0.0;
                float combinedDirY = 0.0;
                float combinedDirZ = 0.0;
                int validPionCount = 0;

                // Sum direction vectors from the two pions (one π+ and one π-)
                if (longPiPlusParticles.size() > 0) {
                    AnaParticleMomB* momPart = dynamic_cast<AnaParticleMomB*>(longPiPlusParticles[0]);
                    if (momPart) {
                        combinedDirX += momPart->DirectionStart[0];
                        combinedDirY += momPart->DirectionStart[1];
                        combinedDirZ += momPart->DirectionStart[2];
                        validPionCount++;
                    }
                }

                if (longPiMinusParticles.size() > 0) {
                    AnaParticleMomB* momPart = dynamic_cast<AnaParticleMomB*>(longPiMinusParticles[0]);
                    if (momPart) {
                        combinedDirX += momPart->DirectionStart[0];
                        combinedDirY += momPart->DirectionStart[1];
                        combinedDirZ += momPart->DirectionStart[2];
                        validPionCount++;
                    }
                }

                // Normalize combined direction vector
                if (validPionCount == 2) { // Only if we have both pions
                    float combinedNorm = sqrt(combinedDirX*combinedDirX + combinedDirY*combinedDirY + combinedDirZ*combinedDirZ);
                    if (combinedNorm > 0) {
                        combinedDirX /= combinedNorm;
                        combinedDirY /= combinedNorm;
                        combinedDirZ /= combinedNorm;
                    }
                }

                // Fill histograms for exactly ONE π+ pion (the first one found)
                if (longPiPlusParticles.size() > 0) {
                    UInt_t piIdx = 0; // Only process the first π+
                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(longPiPlusParticles[piIdx]);
                    if (particlePD) {
                        float piStartX = particlePD->PositionStart[0];
                        float piStartY = particlePD->PositionStart[1];
                        float piStartZ = particlePD->PositionStart[2];

                        if (piStartX > -900 && piStartY > -900 && piStartZ > -900) {
                            // Calculate distance from K0 parent to this π+
                            float piPlusDistance = sqrt((piStartX - parentEndX)*(piStartX - parentEndX) +
                                                       (piStartY - parentEndY)*(piStartY - parentEndY) +
                                                       (piStartZ - parentEndZ)*(piStartZ - parentEndZ));

                            // Fill distance histogram
                            hPiPlusDistance->Fill(piPlusDistance);

                            // Calculate angle between parent->pion vector and pion direction
                            float parentToPionX = piStartX - parentEndX;
                            float parentToPionY = piStartY - parentEndY;
                            float parentToPionZ = piStartZ - parentEndZ;
                            float parentToPionNorm = sqrt(parentToPionX*parentToPionX + parentToPionY*parentToPionY + parentToPionZ*parentToPionZ);

                            if (parentToPionNorm > 0) {
                                // Normalize parent->pion vector
                                parentToPionX /= parentToPionNorm;
                                parentToPionY /= parentToPionNorm;
                                parentToPionZ /= parentToPionNorm;

                                // Get pion direction
                                AnaParticleMomB* momPart = dynamic_cast<AnaParticleMomB*>(particlePD);
                                if (momPart) {
                                    float pionDirX = momPart->DirectionStart[0];
                                    float pionDirY = momPart->DirectionStart[1];
                                    float pionDirZ = momPart->DirectionStart[2];

                                    // Calculate angle between parent->pion and pion direction
                                    float dotProduct = parentToPionX*pionDirX + parentToPionY*pionDirY + parentToPionZ*pionDirZ;
                                    float angle = acos(fabs(dotProduct)) * 180.0 / 3.14159; // Convert to degrees
                                    hPiPlusAngle->Fill(angle);

                                    // Fill 2D histogram: distance vs angle
                                    hPiPlusDistanceVsAngle->Fill(piPlusDistance, angle);

                            // Calculate angle between K0 direction (parent->pion) and pion system direction
                            if (validPionCount == 2) {
                                float combinedDotProduct = parentToPionX*combinedDirX + parentToPionY*combinedDirY + parentToPionZ*combinedDirZ;
                                float combinedAngle = acos(fabs(combinedDotProduct)) * 180.0 / 3.14159; // Convert to degrees
                                hCombinedAngle->Fill(combinedAngle);

                                // Fill 2D histogram: distance vs combined angle
                                hPiPlusDistanceVsCombinedAngle->Fill(piPlusDistance, combinedAngle);
                            }
                                }
                            }
                        }
                    }
                }

                // Fill histograms for exactly ONE π- pion (the first one found)
                if (longPiMinusParticles.size() > 0) {
                    UInt_t piIdx = 0; // Only process the first π-
                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(longPiMinusParticles[piIdx]);
                    if (particlePD) {
                        float piStartX = particlePD->PositionStart[0];
                        float piStartY = particlePD->PositionStart[1];
                        float piStartZ = particlePD->PositionStart[2];

                        if (piStartX > -900 && piStartY > -900 && piStartZ > -900) {
                            // Calculate distance from K0 parent to this π-
                            float piMinusDistance = sqrt((piStartX - parentEndX)*(piStartX - parentEndX) +
                                                        (piStartY - parentEndY)*(piStartY - parentEndY) +
                                                        (piStartZ - parentEndZ)*(piStartZ - parentEndZ));

                            // Fill distance histogram
                            hPiMinusDistance->Fill(piMinusDistance);

                            // Calculate angle between parent->pion vector and pion direction
                            float parentToPionX = piStartX - parentEndX;
                            float parentToPionY = piStartY - parentEndY;
                            float parentToPionZ = piStartZ - parentEndZ;
                            float parentToPionNorm = sqrt(parentToPionX*parentToPionX + parentToPionY*parentToPionY + parentToPionZ*parentToPionZ);

                            if (parentToPionNorm > 0) {
                                // Normalize parent->pion vector
                                parentToPionX /= parentToPionNorm;
                                parentToPionY /= parentToPionNorm;
                                parentToPionZ /= parentToPionNorm;

                                // Get pion direction
                                AnaParticleMomB* momPart = dynamic_cast<AnaParticleMomB*>(particlePD);
                                if (momPart) {
                                    float pionDirX = momPart->DirectionStart[0];
                                    float pionDirY = momPart->DirectionStart[1];
                                    float pionDirZ = momPart->DirectionStart[2];

                                    // Calculate angle between parent->pion and pion direction
                                    float dotProduct = parentToPionX*pionDirX + parentToPionY*pionDirY + parentToPionZ*pionDirZ;
                                    float angle = acos(fabs(dotProduct)) * 180.0 / 3.14159; // Convert to degrees
                                    hPiMinusAngle->Fill(angle);

                                    // Fill 2D histogram: distance vs angle
                                    hPiMinusDistanceVsAngle->Fill(piMinusDistance, angle);

                            // Calculate angle between K0 direction (parent->pion) and pion system direction
                            if (validPionCount == 2) {
                                float combinedDotProduct = parentToPionX*combinedDirX + parentToPionY*combinedDirY + parentToPionZ*combinedDirZ;
                                float combinedAngle = acos(fabs(combinedDotProduct)) * 180.0 / 3.14159; // Convert to degrees

                                // Fill 2D histogram: distance vs combined angle
                                hPiMinusDistanceVsCombinedAngle->Fill(piMinusDistance, combinedAngle);
                            }
                                }
                            }
                        }
                    }
                }
            }

            // End of event processing
            std::cout << "  Event " << i << " processing complete." << std::endl;
            std::cout << std::endl;
        }

        // Progress output every 1000 events
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << " (K^{0}#rightarrow#pi^{+}#pi^{-} events found: " << longPionEventCount << ")" << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total K^{0}#rightarrow#pi^{+}#pi^{-} events found: " << longPionEventCount << std::endl;
    std::cout << "Total π^{+} distances filled: " << hPiPlusDistance->GetEntries() << std::endl;
    std::cout << "Total π^{-} distances filled: " << hPiMinusDistance->GetEntries() << std::endl;
    std::cout << "Total π^{+} angles filled: " << hPiPlusAngle->GetEntries() << std::endl;
    std::cout << "Total π^{-} angles filled: " << hPiMinusAngle->GetEntries() << std::endl;
    std::cout << "Total combined angles filled: " << hCombinedAngle->GetEntries() << std::endl;
    std::cout << std::endl;
    std::cout << "=== Summary ===" << std::endl;
    std::cout << "Found " << longPionEventCount << " events with K0 -> π+ + π- decays" << std::endl;
    std::cout << "All events had pions traveling >10cm and were daughters of K0 based on truth information" << std::endl;

    // Display the distance histograms
    TCanvas *c1 = new TCanvas("c1", "K^{0} Parent to Daughter Distances", 1200, 600);
    c1->Divide(2, 1);

    // π+ distances
    c1->cd(1);
    hPiPlusDistance->Draw("HIST");
    gStyle->SetOptStat(1111);
    gPad->Update();

    // π- distances
    c1->cd(2);
    hPiMinusDistance->Draw("HIST");
    gStyle->SetOptStat(1111);
    gPad->Update();

    c1->Update();
    std::cout << "Distance histograms displayed. Press any key to continue..." << std::endl;
    c1->WaitPrimitive("k");

    // Display the angle histograms
    TCanvas *c2 = new TCanvas("c2", "K^{0} Parent to Daughter Angles", 1200, 600);
    c2->Divide(2, 1);

    // π+ angles
    c2->cd(1);
    hPiPlusAngle->Draw("HIST");
    gStyle->SetOptStat(1111);
    gPad->Update();

    // π- angles
    c2->cd(2);
    hPiMinusAngle->Draw("HIST");
    gStyle->SetOptStat(1111);
    gPad->Update();

    c2->Update();
    std::cout << "Angle histograms displayed. Press any key to continue..." << std::endl;
    c2->WaitPrimitive("k");

    // Display the combined angle histogram
    TCanvas *c3 = new TCanvas("c3", "Combined Angle Distribution", 800, 600);
    hCombinedAngle->Draw("HIST");
    gStyle->SetOptStat(1111);
    gPad->Update();

    c3->Update();
    std::cout << "Combined angle histogram displayed. Press any key to continue..." << std::endl;
    c3->WaitPrimitive("k");

    // Display the 2D histograms: Distance vs Angle
    TCanvas *c4 = new TCanvas("c4", "Distance vs Angle Relationships", 1200, 800);
    c4->Divide(2, 2);

    // π+ distance vs angle
    c4->cd(1);
    hPiPlusDistanceVsAngle->Draw("COLZ");
    gStyle->SetOptStat(1111);
    gPad->Update();

    // π- distance vs angle
    c4->cd(2);
    hPiMinusDistanceVsAngle->Draw("COLZ");
    gStyle->SetOptStat(1111);
    gPad->Update();

    // π+ distance vs combined angle
    c4->cd(3);
    hPiPlusDistanceVsCombinedAngle->Draw("COLZ");
    gStyle->SetOptStat(1111);
    gPad->Update();

    // π- distance vs combined angle
    c4->cd(4);
    hPiMinusDistanceVsCombinedAngle->Draw("COLZ");
    gStyle->SetOptStat(1111);
    gPad->Update();

    c4->Update();
    std::cout << "2D Distance vs Angle histograms displayed. Press any key to continue..." << std::endl;
    c4->WaitPrimitive("k");

    // Clean up
    f->Close();
    delete spill;
}