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

void analyze_k0truerecodist_pion_lengths() {
    // Open the input file
    TFile* file = TFile::Open("../../../../DATA/mc.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Could not open file ../../../../DATA/mc.root" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("ana");
    if (!tree) {
        std::cout << "Error: Could not find tree 'ana' in file" << std::endl;
        return;
    }

    std::cout << "Tree has " << tree->GetEntries() << " entries" << std::endl;

    // Check for required branches
    TBranch* br_k0truerecodist = tree->GetBranch("k0truerecodist");
    TBranch* br_k0dau1recolength = tree->GetBranch("k0dau1recolength");
    TBranch* br_k0dau2recolength = tree->GetBranch("k0dau2recolength");
    TBranch* br_k0dau1truepdg = tree->GetBranch("k0dau1truepdg");
    TBranch* br_k0dau2truepdg = tree->GetBranch("k0dau2truepdg");
    TBranch* br_k0hastrueobject = tree->GetBranch("k0hastrueobject");
    TBranch* br_isk0 = tree->GetBranch("isk0");

    if (!br_k0truerecodist || !br_k0dau1recolength || !br_k0dau2recolength ||
        !br_k0dau1truepdg || !br_k0dau2truepdg) {
        std::cout << "Error: Required branches not found" << std::endl;
        return;
    }

    // Define cumulative length thresholds for pions
    double lengthThresholds[] = {50, 100, 200, 500, 1000};
    int nRanges = sizeof(lengthThresholds)/sizeof(lengthThresholds[0]);

    std::cout << "Creating " << nRanges*nRanges << " histograms for pion length analysis..." << std::endl;

    // Create 2D array of histograms: histograms[i][j] for pion_plus[i], pion_minus[j]
    TH1F*** histograms = new TH1F**[nRanges];
    for (int i = 0; i < nRanges; i++) {
        histograms[i] = new TH1F*[nRanges];
        for (int j = 0; j < nRanges; j++) {
            TString histName = Form("h_k0truerecodist_pion_plus_%.0f_pion_minus_%.0f", lengthThresholds[i], lengthThresholds[j]);
            TString histTitle = Form("k0truerecodist for pion+<%.0fcm, pion-<%.0fcm", lengthThresholds[i], lengthThresholds[j]);
            histograms[i][j] = new TH1F(histName, histTitle, 30, 0, 20);
        }
    }

    // Single pass through data to fill all histograms
    int nProcessed = tree->GetEntries();
    int nValid = 0;
    int nPionPlusFirst = 0;
    int nPionMinusFirst = 0;
    int nNonPion = 0;

    std::cout << "Processing " << nProcessed << " entries..." << std::endl;

    for (Long64_t entry = 0; entry < nProcessed; entry++) {
        tree->GetEntry(entry);

        double k0truerecodist = tree->GetLeaf("k0truerecodist")->GetValue();
        double k0dau1recolength = tree->GetLeaf("k0dau1recolength")->GetValue();
        double k0dau2recolength = tree->GetLeaf("k0dau2recolength")->GetValue();
        int k0dau1truepdg = tree->GetLeaf("k0dau1truepdg")->GetValue();
        int k0dau2truepdg = tree->GetLeaf("k0dau2truepdg")->GetValue();
        int k0hastrueobject = 0;
        int isk0 = 0;

        if (br_k0hastrueobject) {
            k0hastrueobject = tree->GetLeaf("k0hastrueobject")->GetValue();
        }

        if (br_isk0) {
            isk0 = tree->GetLeaf("isk0")->GetValue();
        }

        if (isk0 != 1) continue;

        // Skip invalid entries
        if (k0truerecodist < 0 || k0dau1recolength < 0 || k0dau2recolength < 0) continue;
        // if (br_k0hastrueobject && k0hastrueobject == 0) continue;

        // Identify pion+ (211) and pion- (-211) daughters
        double pion_plus_length = -1;
        double pion_minus_length = -1;

        if (k0dau1truepdg == 211) {
            pion_plus_length = k0dau1recolength;
            pion_minus_length = k0dau2recolength;
            nPionPlusFirst++;
        } else if (k0dau1truepdg == -211) {
            pion_minus_length = k0dau1recolength;
            pion_plus_length = k0dau2recolength;
            nPionMinusFirst++;
        } else if (k0dau2truepdg == 211) {
            pion_plus_length = k0dau2recolength;
            pion_minus_length = k0dau1recolength;
            nPionPlusFirst++;
        } else if (k0dau2truepdg == -211) {
            pion_minus_length = k0dau2recolength;
            pion_plus_length = k0dau1recolength;
            nPionMinusFirst++;
        } else {
            // Neither daughter is a pion
            nNonPion++;
            continue;
        }

        // Skip if we don't have both pion+ and pion-
        if (pion_plus_length < 0 || pion_minus_length < 0) continue;

        // Fill all applicable histograms in one pass
        for (int i = 0; i < nRanges; i++) {
            if (pion_plus_length < lengthThresholds[i]) {
                for (int j = 0; j < nRanges; j++) {
                    if (pion_minus_length < lengthThresholds[j]) {
                        histograms[i][j]->Fill(k0truerecodist);
                    }
                }
            }
        }

        nValid++;
    }

    std::cout << "Processed " << nProcessed << " entries, " << nValid << " valid pion pairs" << std::endl;
    std::cout << "Pion identification statistics:" << std::endl;
    std::cout << "- Pion+ as dau1: " << nPionPlusFirst << " entries" << std::endl;
    std::cout << "- Pion- as dau1: " << nPionMinusFirst << " entries" << std::endl;
    std::cout << "- Non-pion daughters: " << nNonPion << " entries" << std::endl;

    // Create 2D histogram for MPV values
    TH2F* h2d_mpv = new TH2F("h2d_mpv", "MPV of True-Reco Vtx. Distance",
                             nRanges, 0, nRanges, nRanges, 0, nRanges);
    h2d_mpv->GetXaxis()->SetTitle("Pion+ Length Threshold (cm)");
    h2d_mpv->GetYaxis()->SetTitle("Pion- Length Threshold (cm)");
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
                // Try log-normal fit
                TF1* logNormalFit = new TF1(Form("lognormal_pion_%d_%d", i, j), "x > 0 ? [0]*exp(-0.5*((log(x)-[1])/[2])^2)/(x*[2]*sqrt(2*pi)) : 0", 0.01, 20);

                double maxVal = histograms[i][j]->GetMaximum();
                double binMPV = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                double meanVal = histograms[i][j]->GetMean();
                double rmsVal = histograms[i][j]->GetRMS();

                // Better initial parameter estimation for log-normal
                double amplitude = maxVal * 1.2;
                double mu_estimate = log(meanVal);  // mu parameter (mean of log(x))
                double sigma_estimate = TMath::Max(rmsVal/meanVal, 0.1); // sigma parameter (std dev of log(x))

                logNormalFit->SetParameters(amplitude, mu_estimate, sigma_estimate);

                // Parameter limits for log-normal
                logNormalFit->SetParLimits(0, maxVal * 0.1, maxVal * 5.0); // Amplitude
                logNormalFit->SetParLimits(1, -2.0, 2.0);                  // mu (mean of log(x))
                logNormalFit->SetParLimits(2, 0.05, 2.0);                  // sigma (std dev of log(x))

                // Use better fit options for convergence
                int fitResult = histograms[i][j]->Fit(logNormalFit, "QNR", "", 0.01, 20);

                if (fitResult == 0) {
                    // Check fit quality using Chi2/NDF
                    double chi2_ndf = logNormalFit->GetChisquare() / logNormalFit->GetNDF();
                    if (chi2_ndf < 3.0) { // Good fit quality
                        // Extract MPV from log-normal: MPV = exp(mu - sigma^2)
                        double mu = logNormalFit->GetParameter(1);
                        double sigma = logNormalFit->GetParameter(2);
                        mpv_2d = exp(mu - sigma*sigma);
                        // Error propagation for MPV
                        double mu_err = logNormalFit->GetParError(1);
                        double sigma_err = logNormalFit->GetParError(2);
                        mpv_error_2d = mpv_2d * sqrt(mu_err*mu_err + 4*sigma*sigma*sigma_err*sigma_err);
                        nFitted++;
                    } else {
                        // Poor fit quality, use bin MPV instead
                        mpv_2d = binMPV;
                        mpv_error_2d = 0;
                    }
                } else {
                    // Fallback to bin MPV if fit fails
                    mpv_2d = binMPV;
                    mpv_error_2d = 0;
                }

                delete logNormalFit;
            } else if (histograms[i][j]->GetEntries() >= 5) {
                // Use bin MPV for insufficient data but still some entries
                mpv_2d = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                mpv_error_2d = 0;
            } else {
                nInsufficient++;
            }

            // Fill 2D histogram with MPV value (only if MPV <= 2cm)
            if (mpv_2d > 0 && mpv_2d <= 2.0) {
                h2d_mpv->SetBinContent(i+1, j+1, mpv_2d);
            } else if (mpv_2d > 2.0) {
                // Skip entries with MPV > 2cm
                std::cout << "Skipping pion+<" << lengthThresholds[i] << "cm, pion-<" << lengthThresholds[j]
                          << "cm (MPV=" << mpv_2d << "cm > 2cm)" << std::endl;
            }

            // Print progress for significant combinations
            if (histograms[i][j]->GetEntries() >= 10) {
                std::cout << "pion+<" << lengthThresholds[i] << "cm, pion-<" << lengthThresholds[j]
                         << "cm: " << histograms[i][j]->GetEntries() << " entries, MPV=" << mpv_2d << "cm" << std::endl;
            }
        }
    }

    std::cout << "\nLog-Normal fitting summary:" << std::endl;
    std::cout << "- Successfully fitted: " << nFitted << " histograms" << std::endl;
    std::cout << "- Insufficient data: " << nInsufficient << " histograms" << std::endl;

    // Create and save 2D summary plot
    TCanvas* c1 = new TCanvas("c1", "2D MPV vs Pion Lengths", 1000, 800);
    h2d_mpv->Draw("COLZ");

    // Improve appearance
    gStyle->SetPalette(kRainBow);
    h2d_mpv->GetXaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetYaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetZaxis()->SetTitleOffset(1.2);

    c1->SaveAs("mpv_2d_vs_pion_lengths.pdf");
    std::cout << "Saved: mpv_2d_vs_pion_lengths.pdf" << std::endl;

    // Create additional plots for interesting combinations
    std::cout << "\n=== Creating individual plots for all combinations ===" << std::endl;

    // Create plots directory
    gSystem->Exec("mkdir -p plots");

    // Create individual plots for all combinations with significant statistics
    int nSaved = 0;
    for (int i = 0; i < nRanges; i++) {
        for (int j = 0; j < nRanges; j++) {
            int entries = histograms[i][j]->GetEntries();

            // Only save plots with significant statistics
            if (entries >= 20) {
                TCanvas* c_individual = new TCanvas(Form("c_pion_%d_%d", i, j),
                                                   Form("Pion+<%.0fcm, Pion-<%.0fcm", lengthThresholds[i], lengthThresholds[j]),
                                                   800, 600);

                histograms[i][j]->SetLineColor(kBlue);
                histograms[i][j]->SetLineWidth(2);
                histograms[i][j]->SetFillStyle(0);
                histograms[i][j]->Draw("HIST");

                histograms[i][j]->SetTitle(Form("pion+<%.0fcm, pion-<%.0fcm (%d entries)",
                                               lengthThresholds[i], lengthThresholds[j], entries));
                histograms[i][j]->GetXaxis()->SetTitle("k0truerecodist (cm)");
                histograms[i][j]->GetYaxis()->SetTitle("Entries");

                // Add log-normal fit
                TF1* logNormalFit = new TF1(Form("lognormal_pion_%d_%d", i, j), "x > 0 ? [0]*exp(-0.5*((log(x)-[1])/[2])^2)/(x*[2]*sqrt(2*pi)) : 0", 0.01, 20);

                double maxVal = histograms[i][j]->GetMaximum();
                double binMPV = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                double meanVal = histograms[i][j]->GetMean();
                double rmsVal = histograms[i][j]->GetRMS();

                // Better initial parameter estimation for log-normal
                double amplitude = maxVal * 1.2;
                double mu_estimate = log(meanVal);  // mu parameter (mean of log(x))
                double sigma_estimate = TMath::Max(rmsVal/meanVal, 0.1); // sigma parameter (std dev of log(x))

                logNormalFit->SetParameters(amplitude, mu_estimate, sigma_estimate);

                // Parameter limits for log-normal
                logNormalFit->SetParLimits(0, maxVal * 0.1, maxVal * 5.0); // Amplitude
                logNormalFit->SetParLimits(1, -2.0, 2.0);                  // mu (mean of log(x))
                logNormalFit->SetParLimits(2, 0.05, 2.0);                  // sigma (std dev of log(x))

                // Set fit line style BEFORE fitting
                logNormalFit->SetLineColor(kRed);
                logNormalFit->SetLineWidth(3);
                logNormalFit->SetLineStyle(2); // Dashed line

                int fitResult = histograms[i][j]->Fit(logNormalFit, "Q", "", 0.01, 20);

                if (fitResult == 0) {
                    // Check fit quality using Chi2/NDF
                    double chi2_ndf = logNormalFit->GetChisquare() / logNormalFit->GetNDF();
                    double mpv_used = 0;
                    double mpv_error_used = 0;
                    TString mpv_method = "";

                    if (chi2_ndf < 3.0) { // Good fit quality
                        // Extract MPV from log-normal: MPV = exp(mu - sigma^2)
                        double mu = logNormalFit->GetParameter(1);
                        double sigma = logNormalFit->GetParameter(2);
                        mpv_used = exp(mu - sigma*sigma);
                        // Error propagation for MPV
                        double mu_err = logNormalFit->GetParError(1);
                        double sigma_err = logNormalFit->GetParError(2);
                        mpv_error_used = mpv_used * sqrt(mu_err*mu_err + 4*sigma*sigma*sigma_err*sigma_err);
                        mpv_method = "Log-Normal Fit";
                        // Draw the fit line only for good fits
                        logNormalFit->Draw("SAME");
                    } else {
                        // Poor fit quality, use bin MPV instead
                        mpv_used = binMPV;
                        mpv_error_used = 0;
                        mpv_method = "Bin MPV (poor fit)";
                        std::cout << "Poor fit quality (Chi2/NDF=" << chi2_ndf << ") for pion+<" << lengthThresholds[i] << "cm, pion-<" << lengthThresholds[j] << "cm, using bin MPV" << std::endl;
                    }

                    // Add text box with fit parameters - positioned to avoid stats box
                    TPaveText* textBox = new TPaveText(0.15, 0.75, 0.45, 0.9, "NDC");
                    textBox->SetFillColor(0);
                    textBox->SetBorderSize(1);
                    textBox->SetTextSize(0.035);
                    textBox->AddText(Form("MPV = %.3f #pm %.3f", mpv_used, mpv_error_used));
                    textBox->AddText(Form("Entries = %d", entries));
                    textBox->AddText(Form("Method: %s", mpv_method.Data()));
                    textBox->AddText(Form("#chi^{2}/NDF = %.2f", chi2_ndf));
                    textBox->Draw();
                } else {
                    // Fit failed, just draw the histogram
                    std::cout << "Fit failed for pion+<" << lengthThresholds[i] << "cm, pion-<" << lengthThresholds[j] << "cm" << std::endl;
                }

                // Save individual plot
                TString filename = Form("plots/pion_plus_%.0fcm_pion_minus_%.0fcm.pdf", lengthThresholds[i], lengthThresholds[j]);
                c_individual->SaveAs(filename);
                nSaved++;

                delete logNormalFit;
                delete c_individual;
            }
        }
    }

    std::cout << "Saved " << nSaved << " individual plots to plots/ directory" << std::endl;

    // Create summary plot for top combinations
    std::cout << "\n=== Creating summary plot for top combinations ===" << std::endl;

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
        TCanvas* c2 = new TCanvas("c2", "Individual Histograms for Top Pion Combinations", 1200, 800);
        c2->Divide(3, 2);

        for (int p = 0; p < nPlots; p++) {
            int i = entryCounts[p].second.first;
            int j = entryCounts[p].second.second;
            int entries = entryCounts[p].first;

            c2->cd(p + 1);
            histograms[i][j]->SetLineColor(kBlue);
            histograms[i][j]->SetLineWidth(2);
            histograms[i][j]->Draw("HIST");

            histograms[i][j]->SetTitle(Form("pion+<%.0fcm, pion-<%.0fcm (%d entries)",
                                           lengthThresholds[i], lengthThresholds[j], entries));
            histograms[i][j]->GetXaxis()->SetTitle("k0truerecodist (cm)");
            histograms[i][j]->GetYaxis()->SetTitle("Entries");

            // Add log-normal fit if we have enough entries
            if (entries >= 20) {
                TF1* logNormalFit = new TF1(Form("lognormal_pion_top_%d", p), "x > 0 ? [0]*exp(-0.5*((log(x)-[1])/[2])^2)/(x*[2]*sqrt(2*pi)) : 0", 0.01, 20);

                double maxVal = histograms[i][j]->GetMaximum();
                double binMPV = histograms[i][j]->GetBinCenter(histograms[i][j]->GetMaximumBin());
                double meanVal = histograms[i][j]->GetMean();
                double rmsVal = histograms[i][j]->GetRMS();

                // Better initial parameter estimation for log-normal
                double amplitude = maxVal * 1.2;
                double mu_estimate = log(meanVal);  // mu parameter (mean of log(x))
                double sigma_estimate = TMath::Max(rmsVal/meanVal, 0.1); // sigma parameter (std dev of log(x))

                logNormalFit->SetParameters(amplitude, mu_estimate, sigma_estimate);

                // Parameter limits for log-normal
                logNormalFit->SetParLimits(0, maxVal * 0.1, maxVal * 5.0); // Amplitude
                logNormalFit->SetParLimits(1, -2.0, 2.0);                  // mu (mean of log(x))
                logNormalFit->SetParLimits(2, 0.05, 2.0);                  // sigma (std dev of log(x))

                // Set fit line style BEFORE fitting
                logNormalFit->SetLineColor(kRed);
                logNormalFit->SetLineWidth(3);
                logNormalFit->SetLineStyle(2); // Dashed line

                int fitResult = histograms[i][j]->Fit(logNormalFit, "Q", "", 0.01, 20);

                if (fitResult == 0) {
                    // Check fit quality using Chi2/NDF
                    double chi2_ndf = logNormalFit->GetChisquare() / logNormalFit->GetNDF();
                    double mpv_used = 0;
                    double mpv_error_used = 0;
                    TString mpv_method = "";

                    if (chi2_ndf < 1.0) { // Good fit quality
                        // Extract MPV from log-normal: MPV = exp(mu - sigma^2)
                        double mu = logNormalFit->GetParameter(1);
                        double sigma = logNormalFit->GetParameter(2);
                        mpv_used = exp(mu - sigma*sigma);
                        // Error propagation for MPV
                        double mu_err = logNormalFit->GetParError(1);
                        double sigma_err = logNormalFit->GetParError(2);
                        mpv_error_used = mpv_used * sqrt(mu_err*mu_err + 4*sigma*sigma*sigma_err*sigma_err);
                        mpv_method = "Log-Normal Fit";
                        // Draw the fit line only for good fits
                        logNormalFit->Draw("SAME");
                    } else {
                        // Poor fit quality, use bin MPV instead
                        mpv_used = binMPV;
                        mpv_error_used = 0;
                        mpv_method = "Bin MPV (poor fit)";
                    }

                    // Add text box with fit parameters - positioned to avoid stats box
                    TPaveText* textBox = new TPaveText(0.15, 0.75, 0.45, 0.9, "NDC");
                    textBox->SetFillColor(0);
                    textBox->SetBorderSize(1);
                    textBox->SetTextSize(0.035);
                    textBox->AddText(Form("MPV = %.3f #pm %.3f", mpv_used, mpv_error_used));
                    textBox->AddText(Form("Entries = %d", entries));
                    textBox->AddText(Form("Method: %s", mpv_method.Data()));
                    textBox->AddText(Form("#chi^{2}/NDF = %.2f", chi2_ndf));
                    textBox->Draw();
                } else {
                    // Fit failed, just draw the histogram
                    std::cout << "Fit failed for pion+<" << lengthThresholds[i] << "cm, pion-<" << lengthThresholds[j] << "cm" << std::endl;
                }

                delete  logNormalFit;
            }
        }

        c2->SaveAs("top_pion_combinations_individual_histograms.pdf");
        std::cout << "Saved: top_pion_combinations_individual_histograms.pdf" << std::endl;
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
    std::cout << "- mpv_2d_vs_pion_lengths.pdf (2D color plot)" << std::endl;
    std::cout << "- plots/ directory with " << nSaved << " individual 1D histograms" << std::endl;
    if (nPlots > 0) {
        std::cout << "- top_pion_combinations_individual_histograms.pdf (summary of top 6)" << std::endl;
    }
    std::cout << "\nThis analysis identifies pion+ (PDG=211) and pion- (PDG=-211) daughters" << std::endl;
    std::cout << "and creates MPV plots as a function of their individual lengths." << std::endl;
    std::cout << "Individual plots are saved in the plots/ directory with fitted Landau lines." << std::endl;
}
