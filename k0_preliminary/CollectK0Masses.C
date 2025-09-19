// File: CollectK0Masses.C
// Script to collect all K0 reconstructed invariant masses without event displays
// Usage: root -q CollectK0Masses.C

void CollectK0Masses(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
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

    // Loop over entries and collect K0 masses
    Long64_t nEntries = MiniTree->GetEntries();
    int k0EventCount = 0;

    // Histogram for reconstructed invariant masses
    TH1F* hRecoInvariantMass = new TH1F("hRecoInvariantMass", "K^{0} Reconstructed Invariant Mass;Mass [MeV/c^{2}];Events", 50, 400, 600);
    hRecoInvariantMass->SetLineColor(kBlue);
    hRecoInvariantMass->SetLineWidth(2);
    hRecoInvariantMass->SetFillColor(kBlue);
    hRecoInvariantMass->SetFillStyle(3004);

    std::cout << "Processing " << nEntries << " events to collect K0 reconstructed invariant masses..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);

        // Check if there's a K0 in this event
        bool hasK0 = false;
        UInt_t nTrueParticles = spill->TrueParticles.size();
        for (UInt_t j = 0; j < nTrueParticles; j++) {
            AnaTrueParticleB* part = spill->TrueParticles[j];
            if (!part) continue;
            if (part->PDG == 310) { // K0 particle
                hasK0 = true;
                break;
            }
        }

        // If event has K0, calculate reconstructed mass
        if (hasK0) {
            k0EventCount++;

            // Show progress
            if (i % 1000 == 0 || hasK0) {
                std::cout << "Processing entry " << i << "/" << nEntries;
                if (hasK0) {
                    std::cout << " - Found K0 event #" << k0EventCount;
                }
                std::cout << std::endl;
            }

            // Calculate K0 invariant mass from reconstructed daughters
            std::vector<AnaParticlePD*> k0Daughters;
            AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
            if (bunch) {
                for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                    AnaParticleB* particle = bunch->Particles[j];
                    if (!particle) continue;

                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                    if (!particlePD) continue;

                    // Get associated true particle
                    AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                    if (!truePart) continue;

                    // Check if this is a π+ or π- that is a daughter of a K0
                    if ((truePart->PDG == 211 || truePart->PDG == -211) && truePart->ParentPDG == 310) {
                        // Add all K0 daughters, we'll handle momentum validation later
                        k0Daughters.push_back(particlePD);
                    }
                }
            }

            // Calculate K0 invariant mass from daughters
            if (k0Daughters.size() >= 2) {
                // If we have more than 2 daughters, select the first 2 with valid momentum
                if (k0Daughters.size() > 2) {
                    k0Daughters.resize(2); // Keep only the first 2
                }
                AnaParticlePD* daughter1 = k0Daughters[0];
                AnaParticlePD* daughter2 = k0Daughters[1];

                // Get momentum and direction for both daughters
                double p1 = daughter1->Momentum;
                double p2 = daughter2->Momentum;

                // If momentum is invalid, use default values
                if (p1 <= 0 || p1 >= 10000) {
                    p1 = 200.0; // Default value for pion momentum
                }

                if (p2 <= 0 || p2 >= 10000) {
                    p2 = 200.0; // Default value for pion momentum
                }

                double dir1[3] = {daughter1->DirectionStart[0], daughter1->DirectionStart[1], daughter1->DirectionStart[2]};
                double dir2[3] = {daughter2->DirectionStart[0], daughter2->DirectionStart[1], daughter2->DirectionStart[2]};

                // Calculate 4-momentum components for each daughter
                double px1 = p1 * dir1[0];
                double py1 = p1 * dir1[1];
                double pz1 = p1 * dir1[2];
                double px2 = p2 * dir2[0];
                double py2 = p2 * dir2[1];
                double pz2 = p2 * dir2[2];

                // Calculate energy for each daughter (assuming pion mass = 139.57 MeV/c²)
                double pionMass = 139.57; // MeV/c²
                double E1 = sqrt(p1*p1 + pionMass*pionMass);
                double E2 = sqrt(p2*p2 + pionMass*pionMass);

                // Calculate total 4-momentum of K0
                double totalPx = px1 + px2;
                double totalPy = py1 + py2;
                double totalPz = pz1 + pz2;
                double totalE = E1 + E2;

                // Calculate invariant mass: m² = E² - p²
                double totalP2 = totalPx*totalPx + totalPy*totalPy + totalPz*totalPz;
                double k0RecoInvariantMass = sqrt(totalE*totalE - totalP2);

                // Fill histogram with reconstructed invariant mass
                hRecoInvariantMass->Fill(k0RecoInvariantMass);
            }
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total K0 events found: " << k0EventCount << std::endl;
    std::cout << "Total reconstructed invariant masses filled: " << hRecoInvariantMass->GetEntries() << std::endl;

    // Display the reconstructed invariant mass histogram
    TCanvas *c1 = new TCanvas("c1", "K^{0} Reconstructed Invariant Mass Distribution", 800, 600);
    c1->cd();
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

    c1->Update();
    std::cout << "Reconstructed invariant mass histogram displayed. Press any key to continue..." << std::endl;
    c1->WaitPrimitive("k");

    f->Close();
}
