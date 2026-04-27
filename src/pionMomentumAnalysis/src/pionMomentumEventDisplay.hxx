#ifndef pionMomentumEventDisplay_h
#define pionMomentumEventDisplay_h

#include "pdEventDisplay.hxx"

class pionMomentumEventDisplay : public pdEventDisplay {
public:
  pionMomentumEventDisplay();
  virtual ~pionMomentumEventDisplay();

protected:
  virtual void AddAnalysisVariables(OutputManager& output, Int_t tree_index) override;
  virtual void FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) override;
  virtual bool ReadAnalysisData(TTree* tree) override;
  virtual void DrawParticles3D(TEveScene* scene) override;

private:
  static const Int_t kMaxTrjPoints = 50000;

  Int_t _nTrjPoints;
  Int_t _particle_nTrjPoints[kMaxParticles];
  Int_t _particle_trjPointIndex[kMaxParticles];
  Float_t* _trj_X;
  Float_t* _trj_Y;
  Float_t* _trj_Z;
  Float_t* _trj_dEdx;

  enum enumPionMomentumEventDisplayVars {
    edparticle_nTrjPoints = pdEventDisplay::edmaxvars,
    edparticle_trjPointIndex,
    ednTrjPoints,
    edtrj_X,
    edtrj_Y,
    edtrj_Z,
    edtrj_dEdx,
    edpionmaxvars
  };

  ClassDef(pionMomentumEventDisplay, 0);
};

#endif
