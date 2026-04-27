#include "pionMomentumEventDisplay.hxx"

#include "OutputManager.hxx"

#include <TEveElement.h>
#include <TEvePointSet.h>
#include <TEveScene.h>

#include <algorithm>
#include <string>

ClassImp(pionMomentumEventDisplay);

//********************************************************************
pionMomentumEventDisplay::pionMomentumEventDisplay() : pdEventDisplay() {
//********************************************************************
  _eventDisplayClassName = "pionMomentumEventDisplay";
  _nTrjPoints = 0;
  _trj_X = new Float_t[kMaxTrjPoints];
  _trj_Y = new Float_t[kMaxTrjPoints];
  _trj_Z = new Float_t[kMaxTrjPoints];
  _trj_dEdx = new Float_t[kMaxTrjPoints];
}

//********************************************************************
pionMomentumEventDisplay::~pionMomentumEventDisplay() {
//********************************************************************
  delete[] _trj_X;
  delete[] _trj_Y;
  delete[] _trj_Z;
  delete[] _trj_dEdx;
}

//********************************************************************
void pionMomentumEventDisplay::AddAnalysisVariables(OutputManager& output, Int_t tree_index) {
//********************************************************************
  output.AddVectorVar(tree_index, edparticle_nTrjPoints, "ED_particle_nTrjPoints", "I",
                      "Number of trajectory points per particle", ednParticles, "ED_nParticles", -kMaxParticles);
  output.AddVectorVar(tree_index, edparticle_trjPointIndex, "ED_particle_trjPointIndex", "I",
                      "Starting trajectory-point index for each particle", ednParticles, "ED_nParticles", -kMaxParticles);

  output.AddVectorVar(tree_index, edtrj_X, "ED_trj_X", "F", "Trajectory point X positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_Y, "ED_trj_Y", "F", "Trajectory point Y positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_Z, "ED_trj_Z", "F", "Trajectory point Z positions",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
  output.AddVectorVar(tree_index, edtrj_dEdx, "ED_trj_dEdx", "F", "Trajectory point dEdx (placeholder)",
                      ednTrjPoints, "ED_nTrjPoints", -kMaxTrjPoints);
}

//********************************************************************
void pionMomentumEventDisplay::FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
  (void)box;

  Int_t trjIndex = 0;
  Int_t particleIndex = 0;
  for (Int_t i = 0; i < event.nParticles && i < kMaxParticles; ++i) {
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[i]);
    if (!particle) continue;

    const Int_t nTrj = static_cast<Int_t>(particle->TrjPoints.size());
    output.FillVectorVarForceIndex(edparticle_nTrjPoints, nTrj, particleIndex);
    output.FillVectorVarForceIndex(edparticle_trjPointIndex, trjIndex, particleIndex);

    for (Int_t t = 0; t < nTrj && trjIndex < kMaxTrjPoints; ++t) {
      const AnaTrajectoryPointPD& trj = particle->TrjPoints[nTrj - 1 - t];
      output.FillVectorVar(edtrj_X, static_cast<Float_t>(trj.Position.X()));
      output.FillVectorVar(edtrj_Y, static_cast<Float_t>(trj.Position.Y()));
      output.FillVectorVar(edtrj_Z, static_cast<Float_t>(trj.Position.Z()));
      output.FillVectorVar(edtrj_dEdx, static_cast<Float_t>(-999.0));
      output.IncrementCounter(ednTrjPoints);
      ++trjIndex;
    }
    ++particleIndex;
  }
}

//********************************************************************
bool pionMomentumEventDisplay::ReadAnalysisData(TTree* tree) {
//********************************************************************
  tree->SetBranchAddress("ED_nTrjPoints", &_nTrjPoints);
  tree->SetBranchAddress("ED_particle_nTrjPoints", _particle_nTrjPoints);
  tree->SetBranchAddress("ED_particle_trjPointIndex", _particle_trjPointIndex);
  tree->SetBranchAddress("ED_trj_X", _trj_X);
  tree->SetBranchAddress("ED_trj_Y", _trj_Y);
  tree->SetBranchAddress("ED_trj_Z", _trj_Z);
  tree->SetBranchAddress("ED_trj_dEdx", _trj_dEdx);
  return true;
}

//********************************************************************
void pionMomentumEventDisplay::DrawParticles3D(TEveScene* scene) {
//********************************************************************
  pdEventDisplay::DrawParticles3D(scene);
  if (!scene) return;

  TEveElement* particlesElement = nullptr;
  for (TEveElement::List_i it = scene->BeginChildren(); it != scene->EndChildren(); ++it) {
    TEveElement* child = *it;
    if (!child) continue;
    const char* name = child->GetElementName();
    if (name && std::string(name) == "Particles") {
      // Keep the last "Particles" group, which corresponds to the latest draw call.
      particlesElement = child;
    }
  }
  if (!particlesElement) return;

  Int_t trjIndex = 0;
  Int_t particleSlot = 0;
  for (TEveElement::List_i it = particlesElement->BeginChildren(); it != particlesElement->EndChildren(); ++it) {
    if (particleSlot >= _nParticles || particleSlot >= kMaxParticles) break;

    TEveElementList* particleGroup = dynamic_cast<TEveElementList*>(*it);
    if (!particleGroup) {
      ++particleSlot;
      continue;
    }

    const Int_t nTrj = std::max(0, _particle_nTrjPoints[particleSlot]);
    TEveElementList* trjGroup = new TEveElementList("Trajectory Points");
    particleGroup->AddElement(trjGroup);

    if (nTrj > 0) {
      Int_t color = GetParticleColor(_particle_PDG[particleSlot]);
      if (color == kBlack) color = kGray + 1;

      TEvePointSet* firstTrjSet = new TEvePointSet("First TrjPoint");
      firstTrjSet->SetMarkerStyle(29);
      firstTrjSet->SetMarkerSize(2.0);
      firstTrjSet->SetMainColor(color);
      trjGroup->AddElement(firstTrjSet);

      TEvePointSet* trjSet = new TEvePointSet("TrjPoints");
      trjSet->SetMarkerStyle(20);
      trjSet->SetMarkerSize(0.55);
      trjSet->SetMainColor(color);
      trjGroup->AddElement(trjSet);

      for (Int_t t = 0; t < nTrj && trjIndex < _nTrjPoints && trjIndex < kMaxTrjPoints; ++t) {
        const Float_t x = _trj_X[trjIndex];
        const Float_t y = _trj_Y[trjIndex];
        const Float_t z = _trj_Z[trjIndex];
        if (x > -900.f && y > -900.f && z > -900.f) {
          if (t == nTrj - 1) {
            firstTrjSet->SetNextPoint(x, y, z);
          } else {
            trjSet->SetNextPoint(x, y, z);
          }
        }
        ++trjIndex;
      }
    } else {
      trjIndex += nTrj;
    }

    ++particleSlot;
  }
}
