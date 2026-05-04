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
  static const Int_t kMaxMcsSegments = 1000;

  Int_t _nTrjPoints;
  Int_t _particle_nTrjPoints[kMaxParticles];
  Int_t _particle_trjPointIndex[kMaxParticles];
  Float_t* _trj_X;
  Float_t* _trj_Y;
  Float_t* _trj_Z;
  Float_t* _trj_dEdx;

  /// MCS segments computed for the selected MainTrack at fill time.
  /// Chord = straight line between the first/last trajectory points of each segment.
  /// Fit   = PCA-fitted direction line, centered on the segment centroid, length = arc length.
  Int_t _nMcsSegments;
  Int_t _mcs_mainTrack_uniqueID;
  Float_t _mcs_seg_chord_X1[kMaxMcsSegments];
  Float_t _mcs_seg_chord_Y1[kMaxMcsSegments];
  Float_t _mcs_seg_chord_Z1[kMaxMcsSegments];
  Float_t _mcs_seg_chord_X2[kMaxMcsSegments];
  Float_t _mcs_seg_chord_Y2[kMaxMcsSegments];
  Float_t _mcs_seg_chord_Z2[kMaxMcsSegments];
  Float_t _mcs_seg_fit_X1[kMaxMcsSegments];
  Float_t _mcs_seg_fit_Y1[kMaxMcsSegments];
  Float_t _mcs_seg_fit_Z1[kMaxMcsSegments];
  Float_t _mcs_seg_fit_X2[kMaxMcsSegments];
  Float_t _mcs_seg_fit_Y2[kMaxMcsSegments];
  Float_t _mcs_seg_fit_Z2[kMaxMcsSegments];
  /// False when reading EventDisplayData files predating the MCS branches.
  bool _edHasMcsBranches;

  /// MCS from true trajectory (same segmentation as analysis true MCS: SCE-distorted G4 `Position`).
  Int_t _nMcsTrueSegments;
  Int_t _mcs_true_traj_trueID;
  Float_t _mcs_true_seg_chord_X1[kMaxMcsSegments];
  Float_t _mcs_true_seg_chord_Y1[kMaxMcsSegments];
  Float_t _mcs_true_seg_chord_Z1[kMaxMcsSegments];
  Float_t _mcs_true_seg_chord_X2[kMaxMcsSegments];
  Float_t _mcs_true_seg_chord_Y2[kMaxMcsSegments];
  Float_t _mcs_true_seg_chord_Z2[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_X1[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_Y1[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_Z1[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_X2[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_Y2[kMaxMcsSegments];
  Float_t _mcs_true_seg_fit_Z2[kMaxMcsSegments];
  bool _edHasMcsTrueBranches;

  enum enumPionMomentumEventDisplayVars {
    edparticle_nTrjPoints = pdEventDisplay::edmaxvars,
    edparticle_trjPointIndex,
    ednTrjPoints,
    edtrj_X,
    edtrj_Y,
    edtrj_Z,
    edtrj_dEdx,
    ednMcsSegments,
    edmcs_mainTrack_uniqueID,
    edmcs_seg_chord_X1,
    edmcs_seg_chord_Y1,
    edmcs_seg_chord_Z1,
    edmcs_seg_chord_X2,
    edmcs_seg_chord_Y2,
    edmcs_seg_chord_Z2,
    edmcs_seg_fit_X1,
    edmcs_seg_fit_Y1,
    edmcs_seg_fit_Z1,
    edmcs_seg_fit_X2,
    edmcs_seg_fit_Y2,
    edmcs_seg_fit_Z2,
    ednMcsTrueSegments,
    edmcs_true_traj_trueID,
    edmcs_true_seg_chord_X1,
    edmcs_true_seg_chord_Y1,
    edmcs_true_seg_chord_Z1,
    edmcs_true_seg_chord_X2,
    edmcs_true_seg_chord_Y2,
    edmcs_true_seg_chord_Z2,
    edmcs_true_seg_fit_X1,
    edmcs_true_seg_fit_Y1,
    edmcs_true_seg_fit_Z1,
    edmcs_true_seg_fit_X2,
    edmcs_true_seg_fit_Y2,
    edmcs_true_seg_fit_Z2,
    edpionmaxvars
  };

  ClassDef(pionMomentumEventDisplay, 0);
};

#endif
