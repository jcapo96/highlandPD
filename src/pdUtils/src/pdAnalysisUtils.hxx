#ifndef pdAnalysisUtils_h
#define pdAnalysisUtils_h

#include "CutUtils.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "pdDataClasses.hxx"
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
  Float_t ComputePIDA(const AnaParticlePD& part);

  /// Computes the kinetic energy
  Float_t ComputeKineticEnergy(const AnaParticlePD &part);

  /// Compute dedx from dqdx
  Float_t ComputeDeDxFromDqDx(Float_t dqdx, Int_t plane, Float_t x=0, Float_t y=0, Float_t z=0 );

  /// Compute dqdx from dedx
  Float_t ComputeDqDxFromDeDx(Float_t dedx, Int_t plane);

  /// Compute the total electric field at a given position
  Float_t ComputeTotalEField( Float_t x=0, Float_t y=0, Float_t z=0 );
  
  /// Extrapolate the length of a track to a given Z
  Float_t* ExtrapolateToZ(const AnaParticlePD* part, Float_t z, Float_t* posz);

  /// Compute the average dEdx for several resrange bins
  void ComputeBinnedDeDx(const AnaParticlePD* part, Float_t max_resrange, Int_t nbins, Float_t** avg_dedx);

  /// Find the beam true particle 
  AnaTrueParticle* FindBeamTrueParticle(const AnaSpillB& spill);  

  // Add part2 to part1
  void AddParticles(AnaParticlePD* part1, AnaParticlePD* part2);

  // Compute distances between daughters and vertex
  void ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance);
  
  // Get te AnaTrueParticle with a given ID
  AnaTrueParticlePD* GetTrueParticle(AnaEventB* event, Int_t ID);
  AnaTrueParticlePD* GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID);

  // retreieve the BI particle
  AnaParticlePD* GetBeamParticle(const AnaEventC& event);

  // retreieve the true BI particle
  AnaTrueParticlePD* GetTrueBeamParticle(const AnaEventC& event);

  // Fill the counters for several type of true beam daughters
  void FillBeamDaughterCounters(AnaEventB& event, PDCounters& counters);

  //checks if the beam particle selected by Pandora has been correctly selected. Basic implementation for the moment
  bool isBeamLike(AnaParticlePD* part, AnaBeamPD* beam);
  
  // Compute the PID chi2 and ndf for protons
  std::pair< double, int > Chi2PID(const AnaParticlePD& part, TProfile * profile );

  // Methods to compute the beam PDG variables (cannot be used with the piontree since it does not contain TOF and CKOV info)
  std::vector< int > GetPID( const AnaBeamPD& beam, double nominal_momentum );
  PossibleParticleCands2 GetPIDCandidates( const AnaBeamPD& beam, double nominal_momentum );
  PossibleParticleCands2 GetPIDCandidates_CERNCalib( const AnaBeamPD& beam, double nominal_momentum );

  // Compute the track length using the hit positions
  Float_t ComputeTrackLengthFromHitPosition(const AnaParticlePD* part); 

  // Compute the truncated mean of an std vector
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> part); 

  // Calibrated dQdx
  Float_t ComputeCalibrateddQdX(Float_t prim_dqdx, const TVector3& pos);

  // Compute the wire pitch taking into account the direction of the particle
  Float_t Compute3DWirePitch(Int_t planeKey, const TVector3& dir);

}

#endif
