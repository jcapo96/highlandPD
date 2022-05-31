#ifndef CoherentFit_hxx
#define CoherentFit_hxx

#include "pdDataClasses.hxx"
#include "CoherentSample.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"

const int NHITS = 100;
const int NTOYS = 100;

class CoherentFit{
public:
  CoherentFit();
  virtual ~CoherentFit(){}

  CoherentFit(const std::string& filename);

  TTree* GetTreeFromRootFile();
  
  void WriteToRootFile(const std::string& filename);
  
  void CreateCoherentSamples(const Double_t Chi2Cut);
  
  void GenerateTrueMCHistograms(const double RMIN, const double RMAX, const double STEP,
				const double Chi2Cut,
				const double bin_min = 0, const double bin_max = 50, const double bin_width = 0.1,
				const bool normalize = true,
				const bool resize = true);

  void GenerateHistograms(const double RMIN, const double RMAX, const double STEP,
			  const double Chi2Cut,
			  const double bin_min = 0, const double bin_max = 50, const double bin_width = 0.1,
			  const bool normalize = true,
			  const bool resize = true);

  void GenerateTrueBackgroundHistogramsFromTrueSignal(double RMAX = -1);
  
  void GenerateFakeBackground(const double RMIN, const double RMAX, const double STEP,
			      const double Chi2Cut, const double shift_mean, const double shift_sigma);

  void PropagateSystematicErrors();
  void GenerateToySamples();
  void GenerateToyHistograms(CoherentSample* ToySample, const int itoy);
  void InitializeHistogramsForSystematicErrors();
  void FitToySamples();
  
  void SequentialCoherentFit();
  void DataCoherentFit(const CoherentFit* c);

  void SetParametersFromMCFit(const CoherentFit* c);
  
  void NormalizeHistograms();
  void NormalizeHistograms(CoherentSample* cs);

  void SetBackgroundModel(CoherentSample::BackgroundModelEnum m) {fSignalPlusBackground->SetBackgroundModel(m);
                                                                  fBackground->SetBackgroundModel(m);
								  fTrueBackground->SetBackgroundModel(m);}
  
  CoherentSample* GetSignalPlusBackgroundSample() const {return fSignalPlusBackground;}
  CoherentSample* GetSignalSample() const {return fSignal;}
  CoherentSample* GetBackgroundSample() const {return fBackground;}
  CoherentSample* GetTrueSignalSample() const {return fTrueSignal;}
  CoherentSample* GetTrueBackgroundSample() const {return fTrueBackground;}

  void SetSignalPlusBackgroundSample(CoherentSample *s){fSignalPlusBackground = s;}
  void SetSignalSample(CoherentSample *s){fSignal = s;}
  void SetBackgroundSample(CoherentSample *s){fBackground = s;}
  void SetTrueSignalSample(CoherentSample *s){fTrueSignal = s; fSignalPlusBackground->SetTrueSignal(s); fTrueBackground->SetTrueSignal(s);}
  void SetTrueBackgroundSample(CoherentSample *s){fTrueBackground = s;
                                                  fSignalPlusBackground->SetTrueBackground(s); fTrueSignal->SetTrueBackground(s);}

  void ReplaceTrueSignalSample(CoherentSample *s){delete fTrueSignal; SetTrueSignalSample(s);}
  void ReplaceTrueBackgroundSample(CoherentSample *s){delete fTrueBackground; SetTrueBackgroundSample(s);}
  void ReplaceSignalSample(CoherentSample *s){delete fSignal; SetSignalSample(s);}
  void ReplaceBackgroundSample(CoherentSample *s){delete fBackground; SetBackgroundSample(s);}

  bool GetIsMC() {return fIsMC;}
  
private:

  std::string fFilename;
  
  TFile* fFile;
  
  TTree* fTree;
  bool fIsMiniTree;
  bool fIsSystTree;
  bool fIsMC;
  
  AnaSpillB* fSpill;

  CoherentSample* fSignalPlusBackground;  
  CoherentSample* fSignal;
  CoherentSample* fBackground;
  CoherentSample* fTrueSignal;
  CoherentSample* fTrueBackground;
  std::vector<CoherentSample*> fToySamples;

  TH1F* h_toy_A;
  TH1F* h_toy_B;
  TH1F* h_toy_C;
  TH1F* h_toy_D;
  TH1F* h_toy_R;
};

#endif
