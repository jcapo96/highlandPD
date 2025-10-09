// Simple macro to analyze k0truerecodist MPV for different daughter length ranges
// Author: Generated for neutralKaon analysis

#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>

void analyze_k0truerecodist_mpv() {

    // Set ROOT style
    gStyle->SetOptStat(0);

    // Open the input file
    TFile* file = TFile::Open("../../../../DATA/mc.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file ../../../../DATA/mc.root" << std::endl;
        return;
    }

    // Get the analysis tree
    TTree* tree = (TTree*)file->Get("ana");
    if (!tree) {
        std::cerr << "Error: Cannot find tree 'ana' in the file" << std::endl;
        file->Close();
        return;
    }

    std::cout << "Tree has " << tree->GetEntries() << " entries" << std::endl;

    // Check if branches exist first
    TBranch* br_k0truerecodist = tree->GetBranch("k0truerecodist");
    TBranch* br_k0dau1recolength = tree->GetBranch("k0dau1recolength");
    TBranch* br_k0dau2recolength = tree->GetBranch("k0dau2recolength");
    TBranch* br_k0dau1truepdg = tree->GetBranch("k0dau1truepdg");
    TBranch* br_k0dau2truepdg = tree->GetBranch("k0dau2truepdg");
    TBranch* br_isk0 = tree->GetBranch("isk0");
    TBranch* br_k0hastrueobject = tree->GetBranch("k0hastrueobject");

    if (!br_k0truerecodist) {
        std::cerr << "Error: Could not find k0truerecodist branch" << std::endl;
        file->Close();
        return;
    }

    if (!br_k0dau1recolength) {
        std::cerr << "Error: Could not find k0dau1recolength branch" << std::endl;
        file->Close();
        return;
    }

    if (br_isk0) {
        std::cout << "Found isk0 branch - will filter K0 candidates" << std::endl;
    } else {
        std::cout << "Warning: Could not find isk0 branch - processing all entries" << std::endl;
    }

    // Define cumulative length thresholds (each range includes all entries below the threshold)
    double lengthThresholds[] = {20, 30, 40, 60, 80, 100, 120, 150, 200, 300, 400, 500, 750, 1000};
    int nRanges = sizeof(lengthThresholds)/sizeof(lengthThresholds[0]);

    // Create histograms for each cumulative range (0-10cm for k0truerecodist)
    TH1F* histograms[14]; // For dau1 length ranges
    TH1F* histograms_dau2[14]; // For dau2 length ranges
    for (int i = 0; i < nRanges; i++) {
        TString histName = Form("h_k0truerecodist_under_%.0f", lengthThresholds[i]);
        TString histTitle = Form("k0truerecodist for daughter length <%.0f cm", lengthThresholds[i]);
        histograms[i] = new TH1F(histName, histTitle, 20, 0, 10); // 20 bins

        TString histName_dau2 = Form("h_k0truerecodist_dau2_under_%.0f", lengthThresholds[i]);
        TString histTitle_dau2 = Form("k0truerecodist for dau2 length <%.0f cm", lengthThresholds[i]);
        histograms_dau2[i] = new TH1F(histName_dau2, histTitle_dau2, 20, 0, 10); // 20 bins
    }

    // Note: 2D histogram removed as requested

    // Process entries
    int nEntries = tree->GetEntries();
    int nProcessed = 0;
    int nValid = 0;

    std::cout << "Processing " << nEntries << " entries..." << std::endl;

    for (int i = 0; i < nEntries; i++) {
        if (tree->GetEntry(i) <= 0) continue;
        nProcessed++;

        if (nProcessed % 25000 == 0) {
            std::cout << "Processed " << nProcessed << " entries, " << nValid << " valid..." << std::endl;
        }

        // Read data directly from tree leaves (safer approach)
        Float_t k0truerecodist = tree->GetLeaf("k0truerecodist")->GetValue();
        Float_t k0dau1recolength = tree->GetLeaf("k0dau1recolength")->GetValue();
        Float_t k0dau2recolength = tree->GetLeaf("k0dau2recolength")->GetValue();
        Int_t isk0 = -999; // Default value
        Int_t k0hastrueobject = -999; // Default value
        if (br_isk0) {
            isk0 = (Int_t)tree->GetLeaf("isk0")->GetValue();
        }
        if (br_k0hastrueobject) {
            k0hastrueobject = (Int_t)tree->GetLeaf("k0hastrueobject")->GetValue();
        }

        // Filter K0 candidates if isk0 branch exists - try different filtering criteria
        if (br_k0hastrueobject && k0hastrueobject != 1) continue; // Changed from isk0 != 1 to isk0 <= 0

        // Skip invalid values
        if (k0truerecodist < 0 || k0dau1recolength < 0) continue;

        // Note: 2D histogram filling removed as requested

        // Fill histograms for cumulative ranges (each entry goes into all applicable ranges)
        for (int j = 0; j < nRanges; j++) {
            if (k0dau1recolength < lengthThresholds[j]) {
                histograms[j]->Fill(k0truerecodist);
            }
            if (k0dau2recolength < lengthThresholds[j]) {
                histograms_dau2[j]->Fill(k0truerecodist);
            }
        }

        nValid++;
    }

    std::cout << "Processed " << nProcessed << " entries, " << nValid << " valid entries" << std::endl;

    // Colors for different ranges
    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan, kYellow+2, kPink, kGray};

    // Fit histograms and find MPV using Landau fits
    std::cout << "\n=== MPV Analysis Results (Landau Fits) ===" << std::endl;
    std::cout << "Length Range (cm)\tEntries\tBin MPV (cm)\tLandau MPV (cm)" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    bool firstHist = true;
    double maxY = 0;

    std::vector<double> landauMPVs;
    std::vector<double> landauMPVErrors;
    std::vector<double> landauMPVs_dau2;
    std::vector<double> landauMPVErrors_dau2;

    // 2D analysis: MPV vs (dau1_length, dau2_length)
    TH2F* h2d_mpv = new TH2F("h2d_mpv", "MPV of k0truerecodist vs Daughter Lengths",
                             nRanges, 0, nRanges, nRanges, 0, nRanges);
    h2d_mpv->GetXaxis()->SetTitle("DAU1 Length Threshold Index");
    h2d_mpv->GetYaxis()->SetTitle("DAU2 Length Threshold Index");
    h2d_mpv->SetStats(0);

    // Single loop: fit, create canvas, draw, and save
    for (int i = 0; i < nRanges; i++) {
        // Set histogram style
        histograms[i]->SetLineColor(colors[i % 10]);
        histograms[i]->SetLineWidth(2);
        histograms[i]->SetFillStyle(0); // No fill

        // Create individual canvas for this histogram
        TCanvas* c = new TCanvas(Form("c_%d", i), Form("k0truerecodist <%.0f cm", lengthThresholds[i]), 800, 600);
        histograms[i]->SetTitle(Form("k0truerecodist <%.0f cm", lengthThresholds[i]));
        histograms[i]->GetXaxis()->SetTitle("k0truerecodist (cm)");
        histograms[i]->GetYaxis()->SetTitle("Entries");
        histograms[i]->Draw("HIST");

        double landauMPV = 0;
        double landauMPVError = 0;

        if (histograms[i]->GetEntries() < 20) { // Need more entries for reliable fitting
            std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                      << histograms[i]->GetEntries() << "\t\tInsufficient data" << std::endl;
        } else {
            // Find MPV using histogram peak (fallback)
            double binMPV = histograms[i]->GetBinCenter(histograms[i]->GetMaximumBin());

            // Create and fit Landau function
            TF1* landauFit = new TF1(Form("landau_%d", i), "landau", 0, 10);

            // Set initial parameters based on histogram
            double maxVal = histograms[i]->GetMaximum();
            double mpvEstimate = binMPV;
            double sigmaEstimate = (histograms[i]->GetRMS() > 0) ? histograms[i]->GetRMS() : 0.5;

            landauFit->SetParameters(maxVal * 0.8, mpvEstimate, sigmaEstimate);
            landauFit->SetParLimits(0, 0, maxVal * 2);     // Amplitude
            landauFit->SetParLimits(1, 0.5, 5.0);          // MPV: 0.5 to 5.0 cm
            landauFit->SetParLimits(2, 0.1, 2.0);          // Sigma: 0.1 to 2.0 cm

            // Set fit line style BEFORE fitting
            landauFit->SetLineColor(kRed);
            landauFit->SetLineWidth(2);
            landauFit->SetLineStyle(2); // Dashed line

            // Perform fit and draw
            int fitResult = histograms[i]->Fit(landauFit, "QWL", "", 0, 10);

            if (fitResult == 0) { // Successful fit
                landauMPV = landauFit->GetParameter(1); // MPV parameter
                landauMPVError = landauFit->GetParError(1);

                // Explicitly draw the fit line to ensure it appears
                landauFit->Draw("SAME");

                // Create text box with fit parameters
                TPaveText* textBox = new TPaveText(0.65, 0.75, 0.9, 0.9, "NDC");
                textBox->SetFillColor(0);
                textBox->SetBorderSize(1);
                textBox->SetTextSize(0.035);

                double chi2 = landauFit->GetChisquare();
                double ndf = landauFit->GetNDF();

                textBox->AddText(Form("MPV = %.3f #pm %.3f", landauMPV, landauMPVError));
                textBox->AddText(Form("Entries = %d", (int)histograms[i]->GetEntries()));
                if (ndf > 0) {
                    textBox->AddText(Form("#chi^{2}/NDF = %.2f", chi2/ndf));
                }

                textBox->Draw();

                std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                          << histograms[i]->GetEntries() << "\t\t" << binMPV << "\t\t"
                          << landauMPV << " +/- " << landauMPVError << std::endl;
            } else {
                std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                          << histograms[i]->GetEntries() << "\t\t" << binMPV << "\t\tFit Failed" << std::endl;
            }

            delete landauFit;
        }

        landauMPVs.push_back(landauMPV);
        landauMPVErrors.push_back(landauMPVError);

        // Save individual plot
        c->SaveAs(Form("k0truerecodist_under_%.0fcm.pdf", lengthThresholds[i]));
        std::cout << "Saved: k0truerecodist_under_" << lengthThresholds[i] << "cm.pdf" << std::endl;
    }

    std::cout << "\n=== DAU2 ANALYSIS ===" << std::endl;
    std::cout << "Length Range (cm)\tEntries\tBin MPV (cm)\tLandau MPV (cm)" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    // Single loop: fit, create canvas, draw, and save for dau2
    for (int i = 0; i < nRanges; i++) {
        // Set histogram style
        histograms_dau2[i]->SetLineColor(colors[i % 10]);
        histograms_dau2[i]->SetLineWidth(2);
        histograms_dau2[i]->SetFillStyle(0);

        // Calculate bin MPV
        double binMPV = histograms_dau2[i]->GetBinCenter(histograms_dau2[i]->GetMaximumBin());

        // Initialize fit parameters
        double landauMPV = 0;
        double landauMPVError = 0;

        if (histograms_dau2[i]->GetEntries() < 20) {
            std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                      << histograms_dau2[i]->GetEntries() << "\t\t" << binMPV << "\t\tInsufficient data" << std::endl;
        } else {
            // Create Landau fit function
            TF1* landauFit = new TF1(Form("landau_dau2_%d", i), "landau", 0, 10);

            // Estimate initial parameters
            double maxVal = histograms_dau2[i]->GetMaximum();
            double mpvEstimate = binMPV;
            double sigmaEstimate = histograms_dau2[i]->GetRMS();

            // Set initial parameters and limits
            landauFit->SetParameters(maxVal * 0.8, mpvEstimate, sigmaEstimate);
            landauFit->SetParLimits(0, 0, maxVal * 2);     // Amplitude
            landauFit->SetParLimits(1, 0.5, 5.0);          // MPV: 0.5 to 5.0 cm
            landauFit->SetParLimits(2, 0.1, 2.0);          // Sigma: 0.1 to 2.0 cm

            // Set fit line style BEFORE fitting
            landauFit->SetLineColor(kRed);
            landauFit->SetLineWidth(2);
            landauFit->SetLineStyle(2); // Dashed line

            // Perform fit and draw
            int fitResult = histograms_dau2[i]->Fit(landauFit, "QWL", "", 0, 10);

            if (fitResult == 0) { // Successful fit
                landauMPV = landauFit->GetParameter(1); // MPV parameter
                landauMPVError = landauFit->GetParError(1);

                // Explicitly draw the fit line to ensure it appears
                landauFit->Draw("SAME");

                // Create text box with fit parameters
                TPaveText* textBox = new TPaveText(0.65, 0.75, 0.9, 0.9, "NDC");
                textBox->SetFillColor(0);
                textBox->SetBorderSize(1);
                textBox->SetTextSize(0.035);
                textBox->AddText(Form("MPV = %.3f #pm %.3f", landauMPV, landauMPVError));
                textBox->AddText(Form("Entries = %.0f", histograms_dau2[i]->GetEntries()));
                textBox->AddText(Form("#chi^{2}/NDF = %.2f", landauFit->GetChisquare() / landauFit->GetNDF()));
                textBox->Draw();

                std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                          << histograms_dau2[i]->GetEntries() << "\t\t" << binMPV << "\t\t" << landauMPV << std::endl;
            } else {
                std::cout << Form("<%.0f", lengthThresholds[i]) << "\t\t\t"
                          << histograms_dau2[i]->GetEntries() << "\t\t" << binMPV << "\t\tFit Failed" << std::endl;
            }

            delete landauFit;
        }

        landauMPVs_dau2.push_back(landauMPV);
        landauMPVErrors_dau2.push_back(landauMPVError);

        // Create individual canvas for dau2
        TCanvas* c_dau2 = new TCanvas(Form("c_dau2_%d", i), Form("k0truerecodist for dau2 length <%.0f cm", lengthThresholds[i]), 800, 600);
        histograms_dau2[i]->Draw("HIST");
        histograms_dau2[i]->SetTitle(Form("k0truerecodist for dau2 length <%.0f cm", lengthThresholds[i]));
        histograms_dau2[i]->GetXaxis()->SetTitle("k0truerecodist (cm)");
        histograms_dau2[i]->GetYaxis()->SetTitle("Entries");

        // Save individual plot for dau2
        c_dau2->SaveAs(Form("k0truerecodist_dau2_under_%.0fcm.pdf", lengthThresholds[i]));
        std::cout << "Saved: k0truerecodist_dau2_under_" << lengthThresholds[i] << "cm.pdf" << std::endl;
    }

    // Note: 2D plot canvas removed as requested

    // Create summary plot: MPV vs daughter length range
    TCanvas* c3 = new TCanvas("c3", "MPV vs Daughter Length Range (Landau Fits)", 800, 600);

    // Create vectors to store MPV values and range centers
    std::vector<double> mpvValues;
    std::vector<double> rangeCenters;
    std::vector<double> mpvErrors;
    std::vector<double> rangeErrors;

    // Use the Landau MPV values we calculated earlier
    int landauIndex = 0;
    for (int i = 0; i < nRanges; i++) {
        if (histograms[i]->GetEntries() >= 20 && landauIndex < landauMPVs.size()) {
            double mpv = landauMPVs[landauIndex];
            double mpvError = landauMPVErrors[landauIndex];
            double rangeCenter = lengthThresholds[i]; // Center of cumulative range (0 to threshold)
            double rangeError = 0; // No error

            mpvValues.push_back(mpv);
            rangeCenters.push_back(rangeCenter);
            mpvErrors.push_back(mpvError); // Use actual Landau fit error
            rangeErrors.push_back(rangeError);
            landauIndex++;
        }
    }

    // Create TGraphErrors for MPV vs range (with Y-axis error bars)
    TGraphErrors* gr = new TGraphErrors(mpvValues.size(),
                                       rangeCenters.data(),
                                       mpvValues.data(),
                                       nullptr, // No X-axis error bars
                                       mpvErrors.data()); // Y-axis error bars from Landau fits

    gr->SetTitle("MPV of k0truerecodist vs Cumulative Daughter Length Threshold (Landau Fits)");
    gr->GetXaxis()->SetTitle("Daughter Length Threshold (cm)");
    gr->GetYaxis()->SetTitle("MPV of k0truerecodist (cm)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->SetLineColor(kBlue);
    gr->SetMarkerColor(kBlue);
    gr->SetLineWidth(2);

    c3->cd();
    gr->Draw("AP");

    // No linear fit requested

    // Save summary plot
    c3->SaveAs("mpv_vs_daughter_length_range_landau.pdf");

    // Create summary plot for dau2: MPV vs daughter length range
    TCanvas* c4 = new TCanvas("c4", "MPV vs Daughter Length Range (Landau Fits) - DAU2", 800, 600);

    // Create vectors to store MPV values and range centers for dau2
    std::vector<double> mpvValues_dau2;
    std::vector<double> rangeCenters_dau2;
    std::vector<double> mpvErrors_dau2;
    std::vector<double> rangeErrors_dau2;

    for (int i = 0; i < nRanges; i++) {
        if (landauMPVs_dau2[i] > 0) { // Only include successful fits
            mpvValues_dau2.push_back(landauMPVs_dau2[i]);
            rangeCenters_dau2.push_back(lengthThresholds[i] / 2.0);
            mpvErrors_dau2.push_back(landauMPVErrors_dau2[i]);
            rangeErrors_dau2.push_back(lengthThresholds[i] / 2.0);
        }
    }

    // Create TGraphErrors for MPV vs range (with Y-axis error bars) for dau2
    TGraphErrors* gr_dau2 = new TGraphErrors(mpvValues_dau2.size(),
                                           rangeCenters_dau2.data(),
                                           mpvValues_dau2.data(),
                                           nullptr, // No X-axis error bars
                                           mpvErrors_dau2.data()); // Y-axis error bars from Landau fits

    gr_dau2->SetTitle("MPV of k0truerecodist vs Cumulative Daughter Length Threshold (Landau Fits) - DAU2");
    gr_dau2->GetXaxis()->SetTitle("Daughter Length Threshold (cm)");
    gr_dau2->GetYaxis()->SetTitle("MPV of k0truerecodist (cm)");
    gr_dau2->SetMarkerStyle(20);
    gr_dau2->SetMarkerSize(1.5);
    gr_dau2->SetLineColor(kRed);
    gr_dau2->SetMarkerColor(kRed);
    gr_dau2->SetLineWidth(2);

    c4->cd();
    gr_dau2->Draw("AP");

    // Save dau2 summary plot
    c4->SaveAs("mpv_vs_daughter_length_range_landau_dau2.pdf");

    std::cout << "\n=== 2D ANALYSIS: MPV vs (DAU1, DAU2) Length Combinations ===" << std::endl;
    std::cout << "Creating histograms for each dau1-dau2 length combination..." << std::endl;

    // Nested loops: for each dau1 range, iterate over all dau2 ranges
    for (int i = 0; i < nRanges; i++) {  // dau1 index
        for (int j = 0; j < nRanges; j++) {  // dau2 index

            // Create histogram for this specific dau1-dau2 combination
            TString histName = Form("h_k0truerecodist_dau1_%.0f_dau2_%.0f", lengthThresholds[i], lengthThresholds[j]);
            TString histTitle = Form("k0truerecodist for dau1<%.0fcm, dau2<%.0fcm", lengthThresholds[i], lengthThresholds[j]);
            TH1F* hist_2d = new TH1F(histName, histTitle, 20, 0, 10);

            // Fill histogram with data that satisfies both dau1 and dau2 conditions
            for (Long64_t entry = 0; entry < nProcessed; entry++) {
                tree->GetEntry(entry);

                double k0truerecodist = tree->GetLeaf("k0truerecodist")->GetValue();
                double k0dau1recolength = tree->GetLeaf("k0dau1recolength")->GetValue();
                double k0dau2recolength = tree->GetLeaf("k0dau2recolength")->GetValue();
                int k0hastrueobject = tree->GetLeaf("k0hastrueobject")->GetValue();

                if (k0hastrueobject == 0) continue;

                // Check if data satisfies both dau1 and dau2 conditions
                if (k0truerecodist >= 0 && k0dau1recolength >= 0 && k0dau2recolength >= 0 &&
                    k0dau1recolength < lengthThresholds[i] && k0dau2recolength < lengthThresholds[j]) {
                    hist_2d->Fill(k0truerecodist);
                }
            }

            // Calculate MPV for this combination
            double mpv_2d = 0;
            double mpv_error_2d = 0;

            if (hist_2d->GetEntries() >= 20) {
                // Try Landau fit
                TF1* landauFit_2d = new TF1(Form("landau_2d_%d_%d", i, j), "landau", 0, 10);

                double maxVal = hist_2d->GetMaximum();
                double binMPV = hist_2d->GetBinCenter(hist_2d->GetMaximumBin());
                double sigmaEstimate = hist_2d->GetRMS();

                landauFit_2d->SetParameters(maxVal * 0.8, binMPV, sigmaEstimate);
                landauFit_2d->SetParLimits(0, 0, maxVal * 2);
                landauFit_2d->SetParLimits(1, 0.5, 5.0);
                landauFit_2d->SetParLimits(2, 0.1, 2.0);

                int fitResult = hist_2d->Fit(landauFit_2d, "QN", "", 0, 10);

                if (fitResult == 0) {
                    mpv_2d = landauFit_2d->GetParameter(1);
                    mpv_error_2d = landauFit_2d->GetParError(1);
                } else {
                    // Fallback to bin MPV if fit fails
                    mpv_2d = binMPV;
                    mpv_error_2d = 0;
                }

                delete landauFit_2d;
            } else {
                // Use bin MPV for insufficient data
                mpv_2d = hist_2d->GetBinCenter(hist_2d->GetMaximumBin());
                mpv_error_2d = 0;
            }

            // Fill 2D histogram with MPV value
            if (mpv_2d > 0) {
                h2d_mpv->SetBinContent(i+1, j+1, mpv_2d);  // ROOT bins start from 1
            }

            // Set axis labels for 2D histogram
            h2d_mpv->GetXaxis()->SetBinLabel(i+1, Form("%.0f", lengthThresholds[i]));
            h2d_mpv->GetYaxis()->SetBinLabel(j+1, Form("%.0f", lengthThresholds[j]));

            // Print progress for significant combinations
            if (hist_2d->GetEntries() >= 10) {
                std::cout << "dau1<" << lengthThresholds[i] << "cm, dau2<" << lengthThresholds[j]
                         << "cm: " << hist_2d->GetEntries() << " entries, MPV=" << mpv_2d << "cm" << std::endl;
            }

            delete hist_2d;
        }
    }

    // Create and save 2D summary plot
    TCanvas* c5 = new TCanvas("c5", "2D MPV vs Daughter Lengths", 1000, 800);
    h2d_mpv->Draw("COLZ");

    // Add color bar and improve appearance
    gStyle->SetPalette(kRainBow);
    h2d_mpv->GetZaxis()->SetTitle("MPV of k0truerecodist (cm)");
    h2d_mpv->GetXaxis()->SetTitle("DAU1 Length Threshold (cm)");
    h2d_mpv->GetYaxis()->SetTitle("DAU2 Length Threshold (cm)");
    h2d_mpv->GetXaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetYaxis()->SetTitleOffset(1.2);
    h2d_mpv->GetZaxis()->SetTitleOffset(1.2);

    c5->SaveAs("mpv_2d_vs_daughter_lengths.pdf");
    std::cout << "Saved: mpv_2d_vs_daughter_lengths.pdf" << std::endl;

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Plots saved as:" << std::endl;
    std::cout << "\nDAU1 plots:" << std::endl;
    for (int i = 0; i < nRanges; i++) {
        std::cout << "- k0truerecodist_under_" << lengthThresholds[i] << "cm.pdf" << std::endl;
    }
    std::cout << "- mpv_vs_daughter_length_range_landau.pdf" << std::endl;

    std::cout << "\nDAU2 plots:" << std::endl;
    for (int i = 0; i < nRanges; i++) {
        std::cout << "- k0truerecodist_dau2_under_" << lengthThresholds[i] << "cm.pdf" << std::endl;
    }
    std::cout << "- mpv_vs_daughter_length_range_landau_dau2.pdf" << std::endl;

    std::cout << "\n2D Analysis plots:" << std::endl;
    std::cout << "- mpv_2d_vs_daughter_lengths.pdf" << std::endl;

    // Clean up
    file->Close();

    std::cout << "\nAnalysis complete!" << std::endl;
}