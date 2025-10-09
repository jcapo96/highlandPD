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

void k0truerecodist(){

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
    TBranch* br_isk0 = tree->GetBranch("isk0");
    TBranch* br_k0hastrueobject = tree->GetBranch("k0hastrueobject");

    // Create histograms
    TH1F* h_k0truerecodist = new TH1F("h_k0truerecodist", "K0 true recodist", 20, 0, 20);
    h_k0truerecodist->SetXTitle("K0 true recodist");
    h_k0truerecodist->SetYTitle("Events");

    // Fill histograms
    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (i % 25000 == 0) {
            std::cout << "Processing entry " << i << " of " << nentries << std::endl;
            std::cout << i << " is " << br_isk0->GetLeaf("isk0")->GetValue() << std::endl;
            std::cout << i << " has true object " << br_k0hastrueobject->GetLeaf("k0hastrueobject")->GetValue() << std::endl;
        }
        // if (br_isk0->GetLeaf("isk0")->GetValue() != 1) {
        //     continue;
        // }
        // if (br_k0hastrueobject->GetLeaf("k0hastrueobject")->GetValue() == 0) {
        //     continue;
        // }
        h_k0truerecodist->Fill(br_k0truerecodist->GetLeaf("k0truerecodist")->GetValue());
    }

    // Save histograms
    TCanvas* c = new TCanvas("c", "K0 true recodist", 800, 600);
    h_k0truerecodist->Draw();
    c->Update();
    // c->SaveAs("k0truerecodist.png");
    c->WaitPrimitive();
}