#include "dEdxDataCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include "TFile.h"
#include <cassert>

//#define DEBUG

//********************************************************************
dEdxDataCorrection::dEdxDataCorrection(){
//********************************************************************

  std::string filename = std::string(getenv("PROTODUNEEXAMPLEANALYSISROOT"))+"/data/run_5387_Xcalibration.root";

  // read the file with the calibration histograms
  TFile cali_file(filename.c_str());

  // get the calibration histogram
  cali_factor=static_cast<TH1F*>(cali_file.Get("dqdx_X_correction_hist")->Clone());
  cali_factor->SetDirectory(NULL);
}

//********************************************************************
void dEdxDataCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // No correction for data
  if (spill.GetIsMC()) return;

  // Loop over all bunches
  for (unsigned int i = 0; i < spill.Bunches.size(); i++) {
    AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[i]);

    // Loop over all relevant tracks for this variation
    for (UInt_t itrk = 0; itrk<bunch->Particles.size(); itrk++){
      
      AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[itrk]);
      
      // The un-corrected particle
      const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

      if (!original) continue; //?

      for (Int_t i=2;i<3;i++){
        for (UInt_t j=0;j<part->Hits[i].size();j++){
          if (part->Hits[i][j].dQdx<0.1) continue;  // Protection against crazy values
          part->Hits[i][j].dQdx = ComputeCalibratedDqDx(original->Hits[i][j].dQdx,original->Hits[i][j].Position.X());
          part->Hits[i][j].dEdx = pdAnaUtils::ComputeDeDxFromDqDx(part->Hits[i][j].dQdx, i);
          //          std::cout << part->dQdx[plane][j] << " " << original->dQdx[plane][j] << " " << part->dEdx[plane][j] << " " << original->dEdx[plane][j] << " " << part->HitX[plane][j] << std::endl; 
        }
      }
    }
  }
}

//********************************************************************
Float_t dEdxDataCorrection::ComputeCalibratedDqDx(Float_t prim_dqdx, Float_t prim_hitx) {
//********************************************************************
  
  double normalisation_factor=0.983;//for plane 2

  Int_t bin = cali_factor->FindBin(prim_hitx);
  if (bin<1) bin = 1;
  if (bin>148) bin = 148;
  double Cx=cali_factor->GetBinContent(bin);
  double cali_dqdx=Cx*prim_dqdx*normalisation_factor;

  //  std::cout << Cx << " " << prim_hitx << " " << prim_dqdx << " " << cali_dqdx*calib_factor << " " << prim_dedx << " " << cali_dedx  << std::endl; 
  
  return cali_dqdx;
}

