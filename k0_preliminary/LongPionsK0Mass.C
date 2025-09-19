// File: LongPionsK0Mass.C
// Based on original LongPions.C concept: events with long pions (>10cm) + K0 mass calculation
// Only processes events with at least 1 π+ and 1 π- that travel more than 10cm
void LongPionsK0Mass(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
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

    // Loop over entries and collect K0 invariant masses
    Long64_t nEntries = MiniTree->GetEntries();
    int longPionEventCount = 0;

    // Histogram for reconstructed invariant masses
    TH1F* hRecoInvariantMass = new TH1F("hRecoInvariantMass", "K^{0} Reconstructed Invariant Mass;Mass [MeV/c^{2}];Events", 50, 400, 600);
    hRecoInvariantMass->SetLineColor(kBlue);
    hRecoInvariantMass->SetLineWidth(2);
    hRecoInvariantMass->SetFillColor(kBlue);
    hRecoInvariantMass->SetFillStyle(3004);

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);

        // Check if there are long pions (π+ and π- that travel more than 10cm) AND are daughters of K0
        bool hasLongPiPlus = false;
        bool hasLongPiMinus = false;

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

                // Check if this is a π+ or π- daughter of K0 and calculate track length
                if ((truePart->PDG == 211 || truePart->PDG == -211) && truePart->ParentPDG == 310) {
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
                            if (truePart->PDG == 211) {
                                hasLongPiPlus = true;
                            } else if (truePart->PDG == -211) {
                                hasLongPiMinus = true;
                            }
                        }
                    }
                }
            }
        }

        // Only process events with both long π+ and π- daughters of K0
        if (hasLongPiPlus && hasLongPiMinus) {
            longPionEventCount++;

            // Calculate reconstructed K0 invariant mass from daughters
            double k0RecoInvariantMass = -1.0;

            // Find reconstructed π+ and π- daughters of K0
            std::vector<AnaParticleB*> k0Daughters;

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

                    // Check if this is a π+ or π- daughter of K0
                    if ((truePart->PDG == 211 || truePart->PDG == -211) && truePart->ParentPDG == 310) {
                        k0Daughters.push_back(particle);
                    }
                }
            }

            // Calculate invariant mass if we have at least 2 daughters
            if (k0Daughters.size() >= 2) {
                // Use first two daughters (π+ and π-)
                AnaParticleB* daughter1 = k0Daughters[0];
                AnaParticleB* daughter2 = k0Daughters[1];

                // Cast to AnaParticleMomB to access momentum
                AnaParticleMomB* momDaughter1 = dynamic_cast<AnaParticleMomB*>(daughter1);
                AnaParticleMomB* momDaughter2 = dynamic_cast<AnaParticleMomB*>(daughter2);

                if (!momDaughter1 || !momDaughter2) continue;

                // Get momentum for each daughter
                double p1 = momDaughter1->Momentum;
                double p2 = momDaughter2->Momentum;

                // Use default momentum if invalid
                if (p1 <= 0) p1 = 200.0; // MeV/c
                if (p2 <= 0) p2 = 200.0; // MeV/c

                // Get direction for each daughter
                double dir1[3] = {momDaughter1->DirectionStart[0], momDaughter1->DirectionStart[1], momDaughter1->DirectionStart[2]};
                double dir2[3] = {momDaughter2->DirectionStart[0], momDaughter2->DirectionStart[1], momDaughter2->DirectionStart[2]};

                // Calculate 4-momentum components for each daughter
                double pionMass = 139.57; // MeV/c²

                // Daughter 1 (π+)
                double px1 = p1 * dir1[0];
                double py1 = p1 * dir1[1];
                double pz1 = p1 * dir1[2];
                double E1 = sqrt(p1*p1 + pionMass*pionMass);

                // Daughter 2 (π-)
                double px2 = p2 * dir2[0];
                double py2 = p2 * dir2[1];
                double pz2 = p2 * dir2[2];
                double E2 = sqrt(p2*p2 + pionMass*pionMass);

                // Sum 4-momenta to get K0 4-momentum
                double totalPx = px1 + px2;
                double totalPy = py1 + py2;
                double totalPz = pz1 + pz2;
                double totalE = E1 + E2;

                // Calculate invariant mass: m² = E² - p²
                double totalP2 = totalPx*totalPx + totalPy*totalPy + totalPz*totalPz;
                k0RecoInvariantMass = sqrt(totalE*totalE - totalP2);

                // Fill histogram with reconstructed invariant mass
                hRecoInvariantMass->Fill(k0RecoInvariantMass);
            }
        }

        // Progress output every 1000 events
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << " (Long K0 Daughter Pion events found: " << longPionEventCount << ")" << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total long K0 daughter pion events found: " << longPionEventCount << std::endl;
    std::cout << "Total reconstructed invariant masses filled: " << hRecoInvariantMass->GetEntries() << std::endl;

    // Display the reconstructed invariant mass histogram
    TCanvas *c2 = new TCanvas("c2", "K^{0} Reconstructed Invariant Mass Distribution", 800, 600);
    c2->cd();
    hRecoInvariantMass->Draw("HIST");

    // Add statistics box
    gStyle->SetOptStat(1111);
    gPad->Update();

    // Add a line for the expected K0 mass
    TLine *expectedMass = new TLine(497.611, 0, 497.611, hRecoInvariantMass->GetMaximum() * 1.1);
    expectedMass->SetLineColor(kRed);
    expectedMass->SetLineWidth(2);
    expectedMass->SetLineStyle(2);
    expectedMass->Draw();

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hRecoInvariantMass, "Reconstructed Mass", "f");
    leg->AddEntry(expectedMass, "Expected K^{0} Mass", "l");
    leg->Draw();

    c2->Update();
    std::cout << "Reconstructed invariant mass histogram displayed. Press any key to continue..." << std::endl;
    c2->WaitPrimitive("k");

    f->Close();
}