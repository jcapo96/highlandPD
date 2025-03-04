#ifndef CoherentFit_hxx
#define CoherentFit_hxx

#include "pdDataClasses.hxx"
#include "CoherentSample.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom3.h"

const int NHITS = 300;
//const int NTOYS = 100; //not const variable

class CoherentFit{
public:
  CoherentFit();
  virtual ~CoherentFit(){}

  CoherentFit(const std::string& filename, bool IsMC = false);

  TTree* GetTreeFromRootFile();
  TTree* GetSystTree();

  TTree* GetTree() const {return fTree;}
  
  void WriteToRootFile(const std::string& filename);
  
  void CreateCoherentSamples(const Double_t Chi2Cut);
  void CreateSampleLinks();
  
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

  void ComputeSelfSystematicError();
  
  void PropagateSystematicErrors(const std::string& filename, bool apply_toy_weights = true, bool apply_toy_variations = true);
  void ToyLoop(bool apply_toy_weights = true, bool apply_toy_variations = true);
  void GenerateToySample(bool apply_toy_weights = true, bool apply_toy_variations = true, const int itoy = -1);
  void GenerateToyHistograms(CoherentSample* ToySample,
			     bool apply_toy_weights = true, bool apply_toy_variations = true, const int itoy = -1);
  void InitializeHistogramsForSystematicErrors();
  void InitializeHistogramsForSystematicErrors(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> eff);
  void FillSystematicHistograms(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> eff);
  void WriteSystematicHistograms(const std::string& filename);
  void FitToySample(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &eff);
  
  void SequentialCoherentFit(bool minos = true);
  void DataCoherentFit(const CoherentFit* c);
  void ScaleParameters();

  void SetParametersFromMCFit(const CoherentFit* c);
  
  void NormalizeHistograms();
  void NormalizeHistograms(CoherentSample* cs);

  void SetBackgroundModel(CoherentSample::BackgroundModelEnum m) {fSignalPlusBackground->SetBackgroundModel(m);
                                                                  fBackground->SetBackgroundModel(m);
								  fTrueBackground->SetBackgroundModel(m);}
  
  CoherentSample* GetSignalPlusBackgroundSample() const {return fSignalPlusBackground;}
  CoherentSample* GetSignalSample() const {return fSignal;}
  CoherentSample* GetSemiSignalSample() const {return fSemiSignal;}
  CoherentSample* GetBackgroundSample() const {return fBackground;}
  CoherentSample* GetTrueSignalSample() const {return fTrueSignal;}
  CoherentSample* GetTrueSemiSignalSample() const {return fTrueSemiSignal;}
  CoherentSample* GetTrueBackgroundSample() const {return fTrueBackground;}

  CoherentSample* CreateTrueAllKaonsSample(const double RMIN, const double RMAX, const double STEP,
					   const double Chi2Cut,
					   const double bin_min = 0, const double bin_max = 50, const double bin_width = 0.1,
					   const bool normalize = true,
					   const bool resize = true) const;

  CoherentSample* CreateTrueStoppingKaonsSample(const double RMIN, const double RMAX, const double STEP,
						const double Chi2Cut,
						const double bin_min = 0, const double bin_max = 50, const double bin_width = 0.1,
						const bool normalize = true,
						const bool resize = true) const;

  void SetSignalPlusBackgroundSample(CoherentSample *s){fSignalPlusBackground = s;}
  void SetSignalSample(CoherentSample *s){fSignal = s;}
  void SetBackgroundSample(CoherentSample *s){fBackground = s;}
  void SetTrueSignalSample(CoherentSample *s){fTrueSignal = s; fSignalPlusBackground->SetTrueSignal(s); fTrueBackground->SetTrueSignal(s);}
  void SetTrueBackgroundSample(CoherentSample *s){fTrueBackground = s;
                                                  fSignalPlusBackground->SetTrueBackground(s); fTrueSignal->SetTrueBackground(s);}

  void ReplaceTrueSignalSample(CoherentSample *s){delete fTrueSignal; SetTrueSignalSample(s);}
  void ReplaceTrueBackgroundSample(CoherentSample *s){delete fTrueBackground; SetTrueBackgroundSample(s); CreateSampleLinks();}
  void ReplaceSignalSample(CoherentSample *s){delete fSignal; SetSignalSample(s);}
  void ReplaceBackgroundSample(CoherentSample *s){delete fBackground; SetBackgroundSample(s);}

  bool GetIsMC() {return fIsMC;}

  int GetNToys(TTree* t);

  TH1F* GetmpvASystHisto() const {return h_toy_A;}
  TH1F* GetmpvBSystHisto() const {return h_toy_B;}
  TH1F* GetmpvCSystHisto() const {return h_toy_C;}
  
private:

  std::string fFilename;
  
  TFile* fFile;
  
  TTree* fTree;
  TTree* fTreeSystematics;
  bool fIsMiniTree;
  bool fIsSystTree;
  bool fIsMC;
  
  AnaSpillB* fSpill;

  CoherentSample* fSignalPlusBackground;  
  CoherentSample* fSignal;
  CoherentSample* fSemiSignal;
  CoherentSample* fBackground;
  CoherentSample* fSemiBackground;
  CoherentSample* fTrueSignal;
  CoherentSample* fTrueSemiSignal;
  CoherentSample* fTrueBackground;
  CoherentSample* fTrueSemiBackground;
  CoherentSample* fToySample;

  int fNTOYS;
  
  TH1F* h_toy_A;
  TH1F* h_toy_B;
  TH1F* h_toy_C;
  TH1F* h_toy_Eff;
  TH2F* h_toy_dEdx_RR;
};

#endif
