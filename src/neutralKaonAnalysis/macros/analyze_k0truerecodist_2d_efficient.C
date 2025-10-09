#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TBranch.h"
#include <iostream>
#include <vector>

void analyze_k0truerecodist_2d_efficient() {
    // Open the input file
    TFile* file = TFile::Open("../../../../DATA/mc.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Could not open file ../../../../DATA/mc.root" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("ana");
    if (!tree) {
        std::cout << "Error: Could not find tree 'tree' in file" << std::endl;
        return;
    }

    std::cout << "Tree has " << tree->GetEntries() << " entries" << std::endl;

    // Check for required branches
    TBranch* br_k0truerecodist = tree->GetBranch("k0truerecodist");
    TBranch* br_k0dau1recolength = tree->GetBranch("k0dau1recolength");
    TBranch* br_k0dau2recolength = tree->GetBranch("k0dau2recolength");
    TBranch* br_k0hastrueobject = tree->GetBranch("k0hastrueobject");

    if (!br_k0truerecodist || !br_k0dau1recolength || !br_k0dau2recolength) {
        std::cout << "Error: Required branches not found" << std::endl;
        return;
    }

    // Define cumulative length thresholds
    double lengthThresholds[] = {10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200};
    int nRanges = sizeof(lengthThresholds)/sizeof(lengthThresholds[0]);

    std::cout << "Creating " << nRanges*nRanges << " histograms for 2D analysis..." << std::endl;

    // Create 2D array of histograms: histograms[i][j] for dau1[i], dau2[j]
    TH1F*** histograms = new TH1F**[nRanges];
    for (int i = 0; i < nRanges; i++) {
        histograms[i] = new TH1F*[nRanges];
        for (int j = 0; j < nRanges; j++) {
            TString histName = Form("h_k0truerecodist_dau1_%.0f_dau2_%.0f", lengthThresholds[i], lengthThresholds[j]);
            TString histTitle = Form("k0truerecodist for dau1<%.0fcm, dau2<%.0fcm", lengthThresholds[i], lengthThresholds[j]);
            histograms[i][j] = new TH1F(histName, histTitle, 20, 0, 10);
        }
    }

    // Single pass through data to fill all histograms
    int nProcessed = TMath::Min((Long64_t)10000, tree->GetEntries());
    int nValid = 0;

    std::cout << "Processing " << nProcessed << " entries..." << std::endl;

    for (Long64_t entry = 0; entry < nProcessed; entry++) {
        tree->GetEntry(entry);

        double k0truerecodist = tree->GetLeaf("k0truerecodist")->GetValue();
        double k0dau1recolength = tree->GetLeaf("k0dau1recolength")->GetValue();
        double k0dau2recolength = tree->GetLeaf("k0dau2recolength")->GetValue();
        int k0hastrueobject = 0;

        if (br_k0hastrueobject) {
            k0hastrueobject = tree->GetLeaf("k0hastrueobject")->GetValue();
        }

        // Skip invalid entries
        if (k0truerecodist < 0 || k0dau1recolength < 0 || k0dau2recolength < 0) continue;
        // if (br_k0hastrueobject && k0hastrueobject == 0) continue;

        // Fill all applicable histograms in one pass
        for (int i = 0; i < nRanges; i++) {
            if (k0dau1recolength < lengthThresholds[i]) {
                for (int j = 0; j < nRanges; j++) {
                    if (k0dau2recolength < lengthThresholds[j]) {
                        histograms[i][j]->Fill(k0truerecodist);
                    }
                }
            }
        }

        nValid++;
    }

    std::cout << "Processed " << nProcessed << " entries, " << nValid << " valid entries" << std::endl;

    // Create 2D histogram for MPV values
    TH2F* h2d_mpv = new TH2F("h2d_mpv", "MPV of k0truerecodist vs Daughter Lengths",
                             nRanges, 0, nRanges, nRanges, 0, nRanges);
    h2d_mpv->GetXaxis()->SetTitle("DAU1 Length Threshold (cm)");
    h2d_mpv->GetYaxis()->SetTitle("DAU2 Length Threshold (cm)");
    h2d_mpv->GetZaxis()->SetTitle("MPV of k0truerecodist (cm)");
    h2d_mpv->SetStats(0);

    // Set axis labels
    for (int i = 0; i < nRanges; i++) {
        h2d_mpv->GetXaxis()->SetBinLabel(i+1, Form("%.0f", lengthThresholds[i]));
        h2d_mpv->GetYaxis()->SetBinLabel(i+1, Form("%.0f", lengthThresholds[i]));
    }

    std::cout << "\n=== Fitting histograms and extracting MPV values ===" << std::endl;

    // Now fit all histograms and extract MPV values
    int nFitted = 0;
    int nInsufficient = 0;

    for (int i = 0; i < nRanges; i++) {
        for (int j = 0; j < nRanges; j++) {
            double mpv_2d = 0;
            double mpv_error_2d = 0;

            if (histograms[i][j]->GetEntries() >= 20) {
                // Try Landau fit
                TF1* landauFit = new TF1(Form("landau_2d_%d_%d", i, j), "landau", 0, 10);

                double maxVal = histograms[i][j]->GetMaximum();
                double binMPV = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                double sigmaEstimate = histograms[i][j]->GetRMS();

                landauFit->SetParameters(maxVal * 0.8, binMPV, sigmaEstimate);
                landauFit->SetParLimits(0, 0, maxVal * 2);
                landauFit->SetParLimits(1, 0.5, 5.0);
                landauFit->SetParLimits(2, 0.1, 2.0);

                int fitResult = histograms[i][j]->Fit(landauFit, "QN", "", 0, 10);

                if (fitResult == 0) {
                    mpv_2d = landauFit->GetParameter(1);
                    mpv_error_2d = landauFit->GetParError(1);
                    nFitted++;
                } else {
                    // Fallback to bin MPV if fit fails
                    mpv_2d = binMPV;
                    mpv_error_2d = 0;
                }

                delete landauFit;
            } else if (histograms[i][j]->GetEntries() >= 5) {
                // Use bin MPV for insufficient data but still some entries
                mpv_2d = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                mpv_error_2d = 0;
            } else {
                nInsufficient++;
            }

            // Fill 2D histogram with MPV value
            if (mpv_2d > 0) {
                h2d_mpv->SetBinContent(i+1, j+1, mpv_2d);
            }

            // Print progress for significant combinations
            if (histograms[i][j]->GetEntries() >= 10) {
                std::cout << "dau1<" << lengthThresholds[i] << "cm, dau2<" << lengthThresholds[j]
                         << "cm: " << histograms[i][j]->GetEntries() << " entries, MPV=" << mpv_2d << "cm" << std::endl;
            }
        }
    }

    std::cout << "\nFitting summary:" << std::endl;
    std::cout << "- Successfully fitted: " << nFitted << " histograms" << std::endl;
    std::cout << "- Insufficient data: " << nInsufficient << " histograms" << std::endl;

    // Create and save 2D summary plot
    TCanvas* c1 = new TCanvas("c1", "2D MPV vs Daughter Lengths", 1000, 800);
    h2d_mpv->Draw("COLZ");

    // Improve appearance
    gStyle->SetPalette(kRainBow);
    h2d_mpv->GetXaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetYaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetZaxis()->SetTitleOffset(1.2);

    c1->SaveAs("mpv_2d_vs_daughter_lengths_efficient.pdf");
    std::cout << "Saved: mpv_2d_vs_daughter_lengths.pdf" << std::endl;

    // Create additional plots for interesting combinations
    std::cout << "\n=== Creating individual plots for interesting combinations ===" << std::endl;

    // Find combinations with most entries
    std::vector<std::pair<int, std::pair<int, int>>> entryCounts;
    for (int i = 0; i < nRanges; i++) {
        for (int j = 0; j < nRanges; j++) {
            int entries = histograms[i][j]->GetEntries();
            if (entries >= 50) {  // Only plot combinations with significant statistics
                entryCounts.push_back({entries, {i, j}});
            }
        }
    }

    // Sort by entry count (descending)
    std::sort(entryCounts.begin(), entryCounts.end(), std::greater<std::pair<int, std::pair<int, int>>>());

    // Create individual plots for top 6 combinations with fitted lines
    int nPlots = std::min(6, (int)entryCounts.size());
    if (nPlots > 0) {
        TCanvas* c2 = new TCanvas("c2", "Individual Histograms for Top Combinations", 1200, 800);
        c2->Divide(3, 2);

        for (int p = 0; p < nPlots; p++) {
            int i = entryCounts[p].second.first;
            int j = entryCounts[p].second.second;
            int entries = entryCounts[p].first;

            c2->cd(p + 1);
            histograms[i][j]->SetLineColor(kBlue);
            histograms[i][j]->SetLineWidth(2);
            histograms[i][j]->Draw("HIST");

            histograms[i][j]->SetTitle(Form("dau1<%.0fcm, dau2<%.0fcm (%d entries)",
                                           lengthThresholds[i], lengthThresholds[j], entries));
            histograms[i][j]->GetXaxis()->SetTitle("k0truerecodist (cm)");
            histograms[i][j]->GetYaxis()->SetTitle("Entries");

            // Add Landau fit if we have enough entries
            if (entries >= 20) {
                TF1* landauFit = new TF1(Form("landau_top_%d", p), "landau", 0, 10);

                double maxVal = histograms[i][j]->GetMaximum();
                double binMPV = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                double sigmaEstimate = histograms[i][j]->GetRMS();

                landauFit->SetParameters(maxVal * 0.8, binMPV, sigmaEstimate);
                landauFit->SetParLimits(0, 0, maxVal * 2);
                landauFit->SetParLimits(1, 0.5, 5.0);
                landauFit->SetParLimits(2, 0.1, 2.0);

                // Set fit line style
                landauFit->SetLineColor(kRed);
                landauFit->SetLineWidth(2);
                landauFit->SetLineStyle(2); // Dashed line

                int fitResult = histograms[i][j]->Fit(landauFit, "Q", "", 0, 10);

                if (fitResult == 0) {
                    // Fit successful, line is already drawn by the "Q" option
                    // Add text box with fit parameters
                    TPaveText* textBox = new TPaveText(0.65, 0.75, 0.9, 0.9, "NDC");
                    textBox->SetFillColor(0);
                    textBox->SetBorderSize(1);
                    textBox->SetTextSize(0.035);
                    textBox->AddText(Form("MPV = %.3f #pm %.3f", landauFit->GetParameter(1), landauFit->GetParError(1)));
                    textBox->AddText(Form("Entries = %d", entries));
                    textBox->AddText(Form("#chi^{2}/NDF = %.2f", landauFit->GetChisquare() / landauFit->GetNDF()));
                    textBox->Draw();
                } else {
                    // Fit failed, just draw the histogram
                    std::cout << "Fit failed for dau1<" << lengthThresholds[i] << "cm, dau2<" << lengthThresholds[j] << "cm" << std::endl;
                }

                delete landauFit;
            }
        }

        c2->SaveAs("top_combinations_individual_histograms.pdf");
        std::cout << "Saved: top_combinations_individual_histograms.pdf" << std::endl;
    }

    // Clean up
    for (int i = 0; i < nRanges; i++) {
        for (int j = 0; j < nRanges; j++) {
            delete histograms[i][j];
        }
        delete[] histograms[i];
    }
    delete[] histograms;

    file->Close();

    std::cout << "\n=== Analysis Complete ===" << std::endl;
    std::cout << "Generated plots:" << std::endl;
    std::cout << "- mpv_2d_vs_daughter_lengths.pdf" << std::endl;
    if (nPlots > 0) {
        std::cout << "- top_combinations_individual_histograms.pdf" << std::endl;
    }
    std::cout << "\nThis efficient version processes " << nRanges*nRanges << " combinations" << std::endl;
    std::cout << "in a single pass through the data, then fits all histograms." << std::endl;
}
