#ifndef pdAnalysisUtils_h
#define pdAnalysisUtils_h

#include "CutUtils.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "PionAnaDataClasses.hxx"
#include "TProfile.h"

struct PossibleParticleCands2 {

  public:
      bool electron = false;
      bool muon = false;
      bool pion = false;
      bool kaon = false;
      bool proton = false;
      bool deuteron = false;

  enum beamPDG2{
    kElectron = 11,
    kMuon = 13,
    kPion = 211,
    kKaon = 321,
    kProton = 2212,
    kDeuteron = 1000010020
  };



  /*
    inline PossibleParticleCands operator&&(const PossibleParticleCands & b) const {
      return {electron && b.electron, muon && b.muon, pion && b.pion, kaon && b.kaon, proton && b.proton, deuteron && b.deuteron};
    }
    inline PossibleParticleCands operator||(const PossibleParticleCands & b) const {
        return {electron || b.electron, muon || b.muon, pion || b.pion, kaon || b.kaon, proton || b.proton, deuteron || b.deuteron};
    }
    */
    inline operator std::string () const { // overload cast to string
        std::string result = "PossibleParticleCands: [ ";
        if (electron) result += "e ";
        if (muon) result += "mu ";
        if (pion) result += "pi ";
        if (kaon) result += "k ";
        if (proton) result += "p ";
        if (deuteron) result += "e ";
        result += "]";
        return result;
    }
    inline std::vector<int> getPDGCodes() const {
        std::vector<int> result;
        if (electron) result.push_back(kElectron);
        if (muon) result.push_back(kMuon);
        if (pion) result.push_back(kPion);
        if (kaon) result.push_back(kKaon);
        if (proton) result.push_back(kProton);
        if (deuteron) result.push_back(kDeuteron);
        return result;
    }
  };


namespace pdAnaUtils{

  /// Computes the range momentum
  Float_t ComputeRangeMomentum(double trkrange, int pdg);

  /// Computes the CSDARange
  Float_t ComputeCSDARange(double beammom, int pdg);

  /// Compute PIDA
  Float_t ComputePIDA(const AnaParticle& part);

  /// Computes the kinetic energy
  Float_t ComputeKineticEnergy(const AnaParticle &part);

  /// Compute dedx from dqdx
  Float_t ComputeDeDxFromDqDx(Float_t dqdx);

  /// Compute dqdx from dedx
  Float_t ComputeDqDxFromDeDx(Float_t dedx);

  /// Extrapolate the length of a track to a given Z
  Float_t* ExtrapolateToZ(const AnaParticle* part, Float_t z, Float_t* posz);

  /// Compute the average dEdx for several resrange bins
  void ComputeBinnedDeDx(const AnaParticle* part, Float_t max_resrange, Int_t nbins, Float_t** avg_dedx);

  /// Find the beam true particle 
  AnaTrueParticle* FindBeamTrueParticle(const AnaSpillB& spill);  

  // Add part2 to part1
  void AddParticles(AnaParticle* part1, AnaParticle* part2);

  // Compute distances between daughters and vertex
  void ComputeDistanceToVertex(AnaParticle* part, std::vector<Float_t>& distance);
  
  // Get te AnaTrueParticle with a given ID
  AnaTrueParticlePionAna* GetTrueParticle(AnaEventB* event, Int_t ID);
  AnaTrueParticlePionAna* GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID);

  // Fill te counters for several type of true beam daughters
  void FillBeamDaughterCounters(AnaEventB& event, PionAnaCounters& counters);

  //checks if the beam particle selected by Pandora has been correctly selected. Basic implementation for the moment
  bool isBeamLike(AnaParticlePionAna* part, AnaBeamPionAna* beam);
  
  // Compute the PID chi2 and ndf for protons
  std::pair< double, int > Chi2PID(const AnaParticle& part, TProfile * profile );

  // Methods to compute the beam PDG variables (cannot be used with the piontree since it does not contain TOF and CKOV info)
  std::vector< int > GetPID( const AnaBeamPionAna& beam, double nominal_momentum );
  PossibleParticleCands2 GetPIDCandidates( const AnaBeamPionAna& beam, double nominal_momentum );
  PossibleParticleCands2 GetPIDCandidates_CERNCalib( const AnaBeamPionAna& beam, double nominal_momentum );
}

#endif
