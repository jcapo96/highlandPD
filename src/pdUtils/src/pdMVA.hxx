#ifndef pdMVA_h
#define pdMVA_h

#include <map>
#include <vector>
#include <string>

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/Vector3Dfwd.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TVector3.h"
#include "pdDataClasses.hxx"

#include "CalorimetryAlg.hxx"


namespace hl{

  struct MVAPIDResult {
    
    float evalRatio;
    float concentration;
    float coreHaloRatio;
    float conicalness;
    float dEdxStart;
    float dEdxEnd;
    float dEdxEndRatio;
    float length;
    float nSpacePoints;
    unsigned int trackID;
    int isTrack;
    int isStoppingReco;
    
    friend bool           operator <  (const MVAPIDResult & a, const MVAPIDResult & b);
    
    std::map<std::string,double> mvaOutput;
    
  };
  
}

class Track{
public:

  Track(AnaParticleB& part){
    fID = part.UniqueID;
    fEnd = TVector3(part.PositionEnd[0],   part.PositionEnd[1],   part.PositionEnd[2]);
    fLoc = TVector3(part.PositionStart[0], part.PositionStart[1], part.PositionStart[2]);
    fPart = &part;
  }
  
  Int_t fID;
  TVector3 fEnd;
  TVector3 fLoc;
  AnaParticleB* fPart;
  
  Int_t ID(){return fID;}

  template<typename T> inline T Vertex()                        const {return T(fLoc.X(),fLoc.Y(),fLoc.Z()); }
  template<typename T> inline T End()                           const {return T(fEnd.X(),fEnd.Y(),fEnd.Z()); }
  
};

class Shower{
public:
  Int_t fID;
  TVector3 fEnd;
  TVector3 fStart;
  TVector3 fDir;


  Shower(const AnaParticleB& part){
    fID = part.UniqueID;
    fEnd   = TVector3(part.PositionEnd[0],   part.PositionEnd[1],   part.PositionEnd[2]);
    fStart = TVector3(part.PositionStart[0], part.PositionStart[1], part.PositionStart[2]);
    fDir   = TVector3(part.DirectionStart[0], part.DirectionStart[1], part.DirectionStart[2]);
  }

  
  Int_t ID(){return fID;}

  TVector3 ShowerStart(){return fStart;}
  TVector3 Direction(){return   fDir;}
  
};
class SpacePoint{
public:
  Double_t fXYZ[3];

  SpacePoint(const TVector3& pos){
    fXYZ[0] = pos[0];
    fXYZ[1] = pos[1];
    fXYZ[2] = pos[2];
  }
  
  
  const Double_t* XYZ(){return fXYZ;}

};
class Event{};


namespace mvapid {

  //---------------------------------------------------------------
  class MVAAlg {
  public:
    struct SortedObj {
      TVector3 start, end, dir;
      double length;
      std::map<double, AnaHitPD*> hitMap;
    };

    struct SumDistance2 {
      // the TGraph is a data member of the object
      TGraph2D* fGraph;

      SumDistance2(TGraph2D* g) : fGraph(g) {}

      // implementation of the function to be minimized
      double
      operator()(const double* p)
      {

        ROOT::Math::XYZVector x0(p[0], p[2], p[4]);
        ROOT::Math::XYZVector u(p[1], p[3], p[5]);

        u = u.Unit();
        double* x = fGraph->GetX();
        double* y = fGraph->GetY();
        double* z = fGraph->GetZ();
        int npoints = fGraph->GetN();
        double sum = 0;
        for (int i = 0; i < npoints; ++i) {
          ROOT::Math::XYZVector xp(x[i], y[i], z[i]);
          sum += ((xp - x0).Cross(u)).Mag2();
        }
        return sum;
      }
    };

    MVAAlg();

    void Initialize();
    
    void GetDetectorEdges();

    void GetWireNormals();

    void RunPID(AnaEventB& evt,
                std::vector<hl::MVAPIDResult>& result);
                //                art::Assns<recob::Track, hl::MVAPIDResult, void>& trackAssns,
                //                art::Assns<recob::Shower, hl::MVAPIDResult, void>& showerAssns);

  private:
    int IsInActiveVol(const TVector3& pos);

    void PrepareEvent(const AnaEventB& event, const detinfo::DetectorClocksData& clockData);

    void FitAndSortTrack(Track* track, int& isStoppingReco, SortedObj& sortedObj);

    //void SortShower(art::Ptr<recob::Shower> shower,TVector3 dir,int& isStoppingReco,
    //		    mvapid::MVAAlg::SortedObj& sortedShower);
    void SortShower(Shower* shower,
                    int& isStoppingReco,
                    mvapid::MVAAlg::SortedObj& sortedShower);

    void RunPCA(std::vector<AnaHitPD*>& hits,
                std::vector<double>& eVals,
                std::vector<double>& eVecs);

    void _Var_Shape(const SortedObj& track,
                    double& coreHaloRatio,
                    double& concentration,
                    double& conicalness);

    double CalcSegmentdEdxFrac(const detinfo::DetectorClocksData& clock_data,
                               const detinfo::DetectorPropertiesData& det_prop,
                               const SortedObj& track,
                               double start,
                               double end);

    double CalcSegmentdEdxDist(const detinfo::DetectorClocksData& clock_data,
                               const detinfo::DetectorPropertiesData& det_prop,
                               const SortedObj& track,
                               double start,
                               double end);

    double CalcSegmentdEdxDistAtEnd(const detinfo::DetectorClocksData& clock_data,
                                    const detinfo::DetectorPropertiesData& det_prop,
                                    const mvapid::MVAAlg::SortedObj& track,
                                    double distAtEnd);

    int LinFit(Track* track, TVector3& trackPoint, TVector3& trackDir);

    int LinFitShower(Shower* shower,
                     TVector3& showerPoint,
                     TVector3& showerDir);

    const calo::CalorimetryAlg fCaloAlg;

    double fEventT0;

    double fDetMinX, fDetMaxX, fDetMinY, fDetMaxY, fDetMinZ, fDetMaxZ;

    std::map<int, double> fNormToWiresY;
    std::map<int, double> fNormToWiresZ;

    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fHitLabel;
    std::string fSpacePointLabel;
    std::string fTrackingLabel;

    std::vector<Track*> fTracks;
    std::vector<Shower*> fShowers;
    std::vector<SpacePoint*> fSpacePoints;
    std::vector<AnaHitPD*> fHits;

    std::map<Track*, std::vector<AnaHitPD*> > fTracksToHits;
    std::map<Track*, std::vector<SpacePoint*> > fTracksToSpacePoints;
    std::map<Shower*, std::vector<AnaHitPD*> > fShowersToHits;
    std::map<Shower*, std::vector<SpacePoint*> > fShowersToSpacePoints;
    std::map<AnaHitPD*, SpacePoint*> fHitsToSpacePoints;
    //    std::map<SpacePoint*, Hit*> fSpacePointsToHits;

    hl::MVAPIDResult fResHolder;

    TMVA::Reader fReader;

    std::vector<std::string> fMVAMethods;
    std::vector<std::string> fWeightFiles;

    bool fCheatVertex;

    TLorentzVector fVertex4Vect;

  }; // class MVAAlg

} // namespace mvapid

#endif // ifndef pdMVA_h
