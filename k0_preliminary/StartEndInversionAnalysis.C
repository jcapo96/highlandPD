// File: StartEndInversionAnalysis.C
// Analyzes start/end position inversion between true and reconstructed particles
// Counts how many times reconstructed start/end are inverted compared to true particles

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

// Helper function to find reconstructed particle associated with a true particle
AnaParticlePD* FindReconstructedParticle(AnaSpill* spill, Int_t trueParticleID) {
    for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
        AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
        if (!bunch) continue;

        for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
            AnaParticleB* particle = bunch->Particles[j];
            if (!particle) continue;

            // Cast to AnaParticlePD to access ProtoDUNE-specific data
            AnaParticlePD* reconPart = dynamic_cast<AnaParticlePD*>(particle);
            if (!reconPart) continue;

            // Check if this reconstructed particle is associated with the true particle
            AnaTrueParticleB* associatedTrue = reconPart->GetTrueParticle();
            if (associatedTrue && associatedTrue->ID == trueParticleID) {
                return reconPart;
            }
        }
    }
    return nullptr;
}

// Helper function to calculate distance between two 3D points
double CalculateDistance(const Float_t* pos1, const Float_t* pos2) {
    double dx = pos1[0] - pos2[0];
    double dy = pos1[1] - pos2[1];
    double dz = pos1[2] - pos2[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// Helper function to check if start/end positions are inverted
bool IsStartEndInverted(AnaTrueParticleB* truePart, AnaParticlePD* reconPart) {
    if (!truePart || !reconPart) return false;

    // Get true particle positions
    Float_t trueStart[3] = {truePart->Position[0], truePart->Position[1], truePart->Position[2]};
    Float_t trueEnd[3] = {truePart->PositionEnd[0], truePart->PositionEnd[1], truePart->PositionEnd[2]};

    // Get reconstructed particle positions
    Float_t reconStart[3] = {reconPart->PositionStart[0], reconPart->PositionStart[1], reconPart->PositionStart[2]};
    Float_t reconEnd[3] = {reconPart->PositionEnd[0], reconPart->PositionEnd[1], reconPart->PositionEnd[2]};

    // Check if any of the positions are invalid (-999)
    bool trueStartValid = (trueStart[0] != -999 && trueStart[1] != -999 && trueStart[2] != -999);
    bool trueEndValid = (trueEnd[0] != -999 && trueEnd[1] != -999 && trueEnd[2] != -999);
    bool reconStartValid = (reconStart[0] != -999 && reconStart[1] != -999 && reconStart[2] != -999);
    bool reconEndValid = (reconEnd[0] != -999 && reconEnd[1] != -999 && reconEnd[2] != -999);

    // Need both true and reconstructed positions to be valid
    if (!trueStartValid || !trueEndValid || !reconStartValid || !reconEndValid) {
        return false;
    }

    // Calculate distances between true and reconstructed positions
    double distTrueStartToReconStart = CalculateDistance(trueStart, reconStart);
    double distTrueStartToReconEnd = CalculateDistance(trueStart, reconEnd);
    double distTrueEndToReconStart = CalculateDistance(trueEnd, reconStart);
    double distTrueEndToReconEnd = CalculateDistance(trueEnd, reconEnd);

    // Check if reconstructed start is closer to true end than to true start
    // and if reconstructed end is closer to true start than to true end
    bool startInverted = (distTrueStartToReconEnd < distTrueStartToReconStart);
    bool endInverted = (distTrueEndToReconStart < distTrueEndToReconEnd);

    // Both start and end should be inverted for the particle to be considered inverted
    return (startInverted && endInverted);
}

void StartEndInversionAnalysis(const char* dataDir = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA") {

    // Get list of all .root files in the DATA directory
    TSystemDirectory dir(dataDir, dataDir);
    TList *files = dir.GetListOfFiles();

    if (!files) {
        std::cerr << "Error: Could not access directory " << dataDir << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;
    TIter next(files);
    TSystemFile *file;

    while ((file = (TSystemFile*)next())) {
        std::string filename = file->GetName();
        if (filename.length() > 5 && filename.substr(filename.length() - 5) == ".root") {
            // Exclude the reco2 file as requested
            if (filename.find("reco2") != std::string::npos) {
                std::cout << "Skipping excluded file: " << filename << std::endl;
                continue;
            }
            std::string fullPath = std::string(dataDir) + "/" + filename;
            rootFiles.push_back(fullPath);
        }
    }

    if (rootFiles.empty()) {
        std::cerr << "Error: No .root files found in " << dataDir << std::endl;
        return;
    }

    std::cout << "Found " << rootFiles.size() << " root files to process:" << std::endl;
    for (const auto& file : rootFiles) {
        std::cout << "  " << file << std::endl;
    }
    std::cout << std::endl;

    // Global statistics counters
    int totalEvents = 0;
    int totalTrueParticles = 0;
    int totalReconstructedParticles = 0;
    int totalInvertedParticles = 0;

    // Statistics by particle type
    std::map<int, int> trueParticlesByPDG;
    std::map<int, int> reconstructedParticlesByPDG;
    std::map<int, int> invertedParticlesByPDG;

    // Statistics by particle length (true particle length)
    std::map<std::string, int> trueParticlesByLength;
    std::map<std::string, int> reconstructedParticlesByLength;
    std::map<std::string, int> invertedParticlesByLength;

    int fileCounter = 0;
    Long64_t globalEventCounter = 0;

    // Process each file
    for (const auto& filename : rootFiles) {
        fileCounter++;
        std::cout << "Processing file " << fileCounter << "/" << rootFiles.size()
                  << ": " << filename << std::endl;

        // Open the file
        TFile *f = TFile::Open(filename.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file " << filename << std::endl;
            continue;
        }

        // Get the tree
        TTree *MiniTree = (TTree*)f->Get("MiniTree");
        if (!MiniTree) {
            std::cerr << "Could not find MiniTree in file " << filename << std::endl;
            f->Close();
            delete f;
            continue;
        }

        // Set up the object for the branch
        AnaSpill* spill = nullptr;
        MiniTree->SetBranchAddress("Spill", &spill);

        Long64_t nEntries = MiniTree->GetEntries();
        std::cout << "  Events in this file: " << nEntries << std::endl;

        // File-specific counters
        int fileTrueParticles = 0;
        int fileReconstructedParticles = 0;
        int fileInvertedParticles = 0;
        std::map<int, int> fileTrueParticlesByPDG;
        std::map<int, int> fileReconstructedParticlesByPDG;
        std::map<int, int> fileInvertedParticlesByPDG;

        // Process each event
        for (Long64_t i = 0; i < nEntries; i++) {
            globalEventCounter++;
            totalEvents++;

            MiniTree->GetEntry(i);

            if (!spill) continue;

            UInt_t nTrueParticles = spill->TrueParticles.size();

            // Process all true particles
            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Skip particles with invalid positions
                if (truePart->Position[0] == -999 || truePart->Position[1] == -999 || truePart->Position[2] == -999 ||
                    truePart->PositionEnd[0] == -999 || truePart->PositionEnd[1] == -999 || truePart->PositionEnd[2] == -999) {
                    continue;
                }

                totalTrueParticles++;
                fileTrueParticles++;
                trueParticlesByPDG[truePart->PDG]++;
                fileTrueParticlesByPDG[truePart->PDG]++;

                // Calculate true particle length for categorization
                double trueLength = CalculateDistance(truePart->Position, truePart->PositionEnd);
                std::string lengthCategory;
                if (trueLength < 10.0) lengthCategory = "0-10cm";
                else if (trueLength < 50.0) lengthCategory = "10-50cm";
                else if (trueLength < 100.0) lengthCategory = "50-100cm";
                else lengthCategory = "100cm+";

                trueParticlesByLength[lengthCategory]++;

                // Look for associated reconstructed particle
                AnaParticlePD* reconPart = FindReconstructedParticle(spill, truePart->ID);
                if (reconPart) {
                    totalReconstructedParticles++;
                    fileReconstructedParticles++;
                    reconstructedParticlesByPDG[truePart->PDG]++;
                    fileReconstructedParticlesByPDG[truePart->PDG]++;
                    reconstructedParticlesByLength[lengthCategory]++;

                    // Check if start/end positions are inverted
                    if (IsStartEndInverted(truePart, reconPart)) {
                        totalInvertedParticles++;
                        fileInvertedParticles++;
                        invertedParticlesByPDG[truePart->PDG]++;
                        fileInvertedParticlesByPDG[truePart->PDG]++;
                        invertedParticlesByLength[lengthCategory]++;
                    }
                }
            }
        }

        // Print file summary
        std::cout << "  File summary:" << std::endl;
        std::cout << "    True particles: " << fileTrueParticles << std::endl;
        std::cout << "    Reconstructed particles: " << fileReconstructedParticles << std::endl;
        std::cout << "    Inverted particles: " << fileInvertedParticles << std::endl;
        std::cout << "    Inversion rate: " << (fileReconstructedParticles > 0 ? 100.0 * fileInvertedParticles / fileReconstructedParticles : 0.0) << "%" << std::endl;
        std::cout << std::endl;

        f->Close();
        delete f;
    }

    // Print final results
    std::cout << "==========================================" << std::endl;
    std::cout << "START/END INVERSION ANALYSIS SUMMARY" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << "Total true particles: " << totalTrueParticles << std::endl;
    std::cout << "Total reconstructed particles: " << totalReconstructedParticles << std::endl;
    std::cout << "Total inverted particles: " << totalInvertedParticles << std::endl;
    std::cout << "Overall inversion rate: " << (totalReconstructedParticles > 0 ? 100.0 * totalInvertedParticles / totalReconstructedParticles : 0.0) << "%" << std::endl;
    std::cout << std::endl;

    // Print results by particle type
    std::cout << "INVERSION RATE BY PARTICLE TYPE:" << std::endl;
    std::cout << "PDG Code | True | Reco | Inverted | Rate (%)" << std::endl;
    std::cout << "---------|------|------|----------|---------" << std::endl;

    // Get all PDG codes that appear in the data
    std::set<int> allPDGs;
    for (const auto& pair : trueParticlesByPDG) allPDGs.insert(pair.first);
    for (const auto& pair : reconstructedParticlesByPDG) allPDGs.insert(pair.first);

    for (int pdg : allPDGs) {
        int trueCount = trueParticlesByPDG[pdg];
        int recoCount = reconstructedParticlesByPDG[pdg];
        int invertedCount = invertedParticlesByPDG[pdg];
        double rate = (recoCount > 0 ? 100.0 * invertedCount / recoCount : 0.0);

        std::cout << std::setw(8) << pdg << " | "
                  << std::setw(4) << trueCount << " | "
                  << std::setw(4) << recoCount << " | "
                  << std::setw(8) << invertedCount << " | "
                  << std::setw(8) << std::fixed << std::setprecision(2) << rate << std::endl;
    }
    std::cout << std::endl;

    // Print results by particle length
    std::cout << "INVERSION RATE BY PARTICLE LENGTH:" << std::endl;
    std::cout << "Length Range | True | Reco | Inverted | Rate (%)" << std::endl;
    std::cout << "-------------|------|------|----------|---------" << std::endl;

    std::vector<std::string> lengthCategories = {"0-10cm", "10-50cm", "50-100cm", "100cm+"};
    for (const std::string& category : lengthCategories) {
        int trueCount = trueParticlesByLength[category];
        int recoCount = reconstructedParticlesByLength[category];
        int invertedCount = invertedParticlesByLength[category];
        double rate = (recoCount > 0 ? 100.0 * invertedCount / recoCount : 0.0);

        std::cout << std::setw(12) << category << " | "
                  << std::setw(4) << trueCount << " | "
                  << std::setw(4) << recoCount << " | "
                  << std::setw(8) << invertedCount << " | "
                  << std::setw(8) << std::fixed << std::setprecision(2) << rate << std::endl;
    }
    std::cout << std::endl;

    // Print detailed results for common particle types
    std::cout << "DETAILED RESULTS FOR COMMON PARTICLES:" << std::endl;
    std::vector<int> commonPDGs = {211, -211, 111, 13, -13, 11, -11, 321, -321, 2212, -2212};
    for (int pdg : commonPDGs) {
        if (trueParticlesByPDG.find(pdg) != trueParticlesByPDG.end()) {
            int trueCount = trueParticlesByPDG[pdg];
            int recoCount = reconstructedParticlesByPDG[pdg];
            int invertedCount = invertedParticlesByPDG[pdg];
            double rate = (recoCount > 0 ? 100.0 * invertedCount / recoCount : 0.0);

            std::string particleName;
            switch (pdg) {
                case 211: particleName = "π+"; break;
                case -211: particleName = "π-"; break;
                case 111: particleName = "π0"; break;
                case 13: particleName = "μ+"; break;
                case -13: particleName = "μ-"; break;
                case 11: particleName = "e+"; break;
                case -11: particleName = "e-"; break;
                case 321: particleName = "K+"; break;
                case -321: particleName = "K-"; break;
                case 2212: particleName = "p"; break;
                case -2212: particleName = "p̄"; break;
                default: particleName = "PDG" + std::to_string(pdg); break;
            }

            std::cout << "  " << particleName << " (" << pdg << "): "
                      << invertedCount << "/" << recoCount << " ("
                      << std::fixed << std::setprecision(2) << rate << "%)" << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "Analysis completed successfully!" << std::endl;
}
