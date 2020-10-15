#include "pdMVA.hxx"

#include "TPrincipal.h"
#include <Fit/Fitter.h>
#include "Math/Functor.h"
#include "AnalysisUtils.hxx"

#include <cmath>


bool operator < (const hl::MVAPIDResult & a, const hl::MVAPIDResult & b)
  {
    return a.nSpacePoints<b.nSpacePoints;
  }


//**********************************************************
//mvapid::MVAAlg::MVAAlg(): fCaloAlg("CalorimetryAlg"),fReader(""){
mvapid::MVAAlg::MVAAlg(): fReader(""){
//**********************************************************

  /*
  fHitLabel = pset.get<std::string>("HitLabel");
  fTrackLabel = pset.get<std::string>("TrackLabel");
  fShowerLabel = pset.get<std::string>("ShowerLabel");
  fSpacePointLabel = pset.get<std::string>("SpacePointLabel");
  fTrackingLabel = pset.get<std::string>("TrackingLabel", "");

  fCheatVertex = pset.get<bool>("CheatVertex", false);
  */
  fReader.AddVariable("evalRatio",    &fResHolder.evalRatio);
  fReader.AddVariable("coreHaloRatio",&fResHolder.coreHaloRatio);
  fReader.AddVariable("concentration",&fResHolder.concentration);
  fReader.AddVariable("conicalness",  &fResHolder.conicalness);
  fReader.AddVariable("dEdxStart",    &fResHolder.dEdxStart);
  fReader.AddVariable("dEdxEnd",      &fResHolder.dEdxEnd);
  fReader.AddVariable("dEdxEndRatio", &fResHolder.dEdxEndRatio);

  // TODO
  //  fMVAMethods.push_back("ANN");
  /*
  fMVAMethods.push_back("electron");
  fMVAMethods.push_back("muon");
  fMVAMethods.push_back("proton");
  fMVAMethods.push_back("photon");
  fMVAMethods.push_back("pich");

  //  fMVAMethods = pset.get<std::vector<std::string>>("MVAMethods");
  //  std::vector<std::string> weightFileBnames;// = pset.get<std::vector<std::string>>("WeightFiles");
  //  fWeightFiles.push_back("electron_all_ANN.weights.xml");
  fWeightFiles.push_back("electron_all_BDT.weights.xml");
  fWeightFiles.push_back("muon_all_BDT.weights.xml");
  fWeightFiles.push_back("proton_all_BDT.weights.xml");
  fWeightFiles.push_back("photon_all_BDT.weights.xml");
  fWeightFiles.push_back("pich_all_BDT.weights.xml");
  */

  fMVAMethods.push_back("electron");
  fWeightFiles.push_back("cnn_emtrkmichl_pitch_5_wire_48_drift_48_down_6_mean_notes_protoneBeamAndCosmicsMCC11.pb");

  
  //  weightFileBnames.push_back("weights/MuEMVA_ANN.weights.xml");
  //  weightFileBnames.push_back("weights/MuEMVA_BDT.weights.xml");
  //    fWeightFiles.push_back("weights/MuEMVA_ANN.weights.xml");
  //  fWeightFiles.push_back("weights/MuEMVA_BDT.weights.xml");
  
  if (fMVAMethods.size() != fWeightFiles.size()) {
    std::cerr << "Mismatch in number of MVA methods and weight files!" << std::endl;
    exit(1);
  }

  for (unsigned int iMethod = 0; iMethod != fMVAMethods.size(); ++iMethod) {
    fReader.BookMVA(fMVAMethods[iMethod], fWeightFiles[iMethod]);
  }
}


//**********************************************************
int mvapid::MVAAlg::IsInActiveVol(const TVector3& pos){
//**********************************************************
  
  const double fiducialDist = 5.0;

  if (pos.X() > (fDetMinX + fiducialDist) && pos.X() < (fDetMaxX - fiducialDist) &&
      pos.Y() > (fDetMinY + fiducialDist) && pos.Y() < (fDetMaxY - fiducialDist) &&
      pos.Z() > (fDetMinZ + fiducialDist) && pos.Z() < (fDetMaxZ - fiducialDist))
    return 1;
  else
    return 0;
}

//**********************************************************
void mvapid::MVAAlg::GetDetectorEdges(){
//**********************************************************
 
  // Taken from PionAnalyzer_module output;
  fDetMinX = -376.85 ;
  fDetMaxX = 376.85  ;
  fDetMinY = 0       ;
  fDetMaxY = 607.499 ;
  fDetMinZ = -0.49375;
  fDetMaxZ = 695.286 ;

  /*
  //  art::ServiceHandle<geo::Geometry const> geom;

  for (unsigned int t = 0; t < geom->TotalNTPC(); t++) {
    if (geom->TPC(t).MinX() < fDetMinX) fDetMinX = geom->TPC(t).MinX();
    if (geom->TPC(t).MaxX() > fDetMaxX) fDetMaxX = geom->TPC(t).MaxX();
    if (geom->TPC(t).MinY() < fDetMinY) fDetMinY = geom->TPC(t).MinY();
    if (geom->TPC(t).MaxY() > fDetMaxY) fDetMaxY = geom->TPC(t).MaxY();
    if (geom->TPC(t).MinZ() < fDetMinZ) fDetMinZ = geom->TPC(t).MinZ();
    if (geom->TPC(t).MaxZ() > fDetMaxZ) fDetMaxZ = geom->TPC(t).MaxZ();
  }
  */
}

//**********************************************************
void mvapid::MVAAlg::GetWireNormals(){
//**********************************************************  

  //Get normals to wires for each plane in the detector
  //This assumes equal numbers of TPCs in each cryostat and equal numbers of planes in each TPC
  
  fNormToWiresY.clear();
  fNormToWiresZ.clear();

  int planeKey;

  // Numbers from PionAnalyzer_module output
  double dirY_0 = 0.812012;
  double dirZ_0 = 0.58364;

  int NTPC=12;
  planeKey=0;
  for (int t=0;t<NTPC;t++){
    for (int p=0;p<3;p++){

      double dirY=0;
      double dirZ=0;

      if (p==0){
        if (t%2 == 0) dirZ = -dirZ_0;
        else          dirZ =  dirZ_0;
        dirY = dirY_0;
      }
      else if (p==1){
        if (t%2 == 0) dirZ = -dirZ_0;
        else          dirZ =  dirZ_0;
        dirY = -dirY_0;
      }
      else if (p==2){
        if (t%2 == 0) dirY =  -1;
        else          dirY =   1;
        dirZ = 0;
      }
       
      fNormToWiresY.insert(std::make_pair(planeKey, -dirZ)); //y component of normal
      fNormToWiresZ.insert(std::make_pair(planeKey,  dirY)); //z component of normal

      planeKey++;
    }
  }

  
  /*
    art::ServiceHandle<geo::Geometry const> geom;

  for (geo::PlaneGeo const& plane : geom->IteratePlanes()) {
    std::string id = std::string(plane.ID());
    int pcryo = id.find("C");
    int ptpc = id.find("T");
    int pplane = id.find("P");
    std::string scryo = id.substr(pcryo + 2, 2);
    std::string stpc = id.substr(ptpc + 2, 2);
    std::string splane = id.substr(pplane + 2, 2);
    int icryo = std::stoi(scryo);
    int itpc = std::stoi(stpc);
    int iplane = std::stoi(splane);
    planeKey = icryo * geom->NTPC(0) * geom->Nplanes(0, 0) + itpc * geom->Nplanes(0, 0) +
               iplane; //single index for all planes in detector
    fNormToWiresY.insert(
      std::make_pair(planeKey, -plane.Wire(0).Direction().Z())); //y component of normal
    fNormToWiresZ.insert(
      std::make_pair(planeKey, plane.Wire(0).Direction().Y())); //z component of normal
  }
  */


}

//**********************************************************
void mvapid::MVAAlg::Initialize(){
//**********************************************************

  GetDetectorEdges();
  GetWireNormals();
}

//**********************************************************
void mvapid::MVAAlg::RunPID(AnaEventB& evt, std::vector<hl::MVAPIDResult>& result){
//**********************************************************
  //                       art::Assns<recob::Track, hl::MVAPIDResult, void>& trackAssns,
  //                       art::Assns<recob::Shower, hl::MVAPIDResult, void>& showerAssns)


  
  detinfo::DetectorClocksData clockData;
  detinfo::DetectorPropertiesData detProp;
  /*
  auto const clockData =
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  */


  
  this->PrepareEvent(evt, clockData);

  
  for (auto trackIter = fTracks.begin(); trackIter != fTracks.end(); ++trackIter) {
    mvapid::MVAAlg::SortedObj sortedObj;

    std::vector<double> eVals, eVecs;
    int isStoppingReco;

    this->RunPCA(fTracksToHits[*trackIter], eVals, eVecs);

    double evalRatio;
    if (eVals[0] < 0.0001)
      evalRatio = 0.0;
    else
      evalRatio = std::sqrt(eVals[1] * eVals[1] + eVals[2] * eVals[2]) / eVals[0];
    this->FitAndSortTrack(*trackIter, isStoppingReco, sortedObj);
    double coreHaloRatio, concentration, conicalness;
    this->_Var_Shape(sortedObj, coreHaloRatio, concentration, conicalness);
    double dEdxStart = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0., 0.05);
    double dEdxEnd = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0.9, 1.0);
    double dEdxPenultimate = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0.8, 0.9);

    fResHolder.isTrack = 1;
    fResHolder.isStoppingReco = isStoppingReco;
    fResHolder.nSpacePoints = sortedObj.hitMap.size();
    fResHolder.trackID = (*trackIter)->ID();
    fResHolder.evalRatio = evalRatio;
    fResHolder.concentration = concentration;
    fResHolder.coreHaloRatio = coreHaloRatio;
    fResHolder.conicalness = conicalness;
    fResHolder.dEdxStart = dEdxStart;
    fResHolder.dEdxEnd = dEdxEnd;
    if (dEdxPenultimate < 0.1)
      fResHolder.dEdxEndRatio = 1.0;
    else
      fResHolder.dEdxEndRatio = dEdxEnd / dEdxPenultimate;
    fResHolder.length = sortedObj.length;

    std::cout << "---------------------" << std::endl;
    static_cast<AnaParticlePD*>((*trackIter)->fPart)->Print();
    for (auto methodIter = fMVAMethods.begin(); methodIter != fMVAMethods.end(); ++methodIter) {
      fResHolder.mvaOutput[*methodIter] = fReader.EvaluateMVA(*methodIter);
      std::cout << "MVA output " << *methodIter << " " << fResHolder.mvaOutput[*methodIter] << " " << static_cast<AnaParticlePD*>((*trackIter)->fPart)->CNNscore[0]<< std::endl;
      
    }
    result.push_back(fResHolder);
    //    util::CreateAssn(evt, result, *trackIter, trackAssns);
  }

  for (auto showerIter = fShowers.begin(); showerIter != fShowers.end(); ++showerIter) {
    mvapid::MVAAlg::SortedObj sortedObj;

    std::vector<double> eVals, eVecs;
    int isStoppingReco;

    this->RunPCA(fShowersToHits[*showerIter], eVals, eVecs);

    double evalRatio;
    if (eVals[0] < 0.0001)
      evalRatio = 0.0;
    else
      evalRatio = std::sqrt(eVals[1] * eVals[1] + eVals[2] * eVals[2]) / eVals[0];

    this->SortShower(*showerIter, isStoppingReco, sortedObj);

    double coreHaloRatio, concentration, conicalness;
    this->_Var_Shape(sortedObj, coreHaloRatio, concentration, conicalness);
    double dEdxStart = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0., 0.05);
    double dEdxEnd = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0.9, 1.0);
    double dEdxPenultimate = CalcSegmentdEdxFrac(clockData, detProp, sortedObj, 0.8, 0.9);

    fResHolder.isTrack = 0;
    fResHolder.isStoppingReco = isStoppingReco;
    fResHolder.nSpacePoints = sortedObj.hitMap.size();
    fResHolder.trackID =
      (*showerIter)->ID() + 1000; //For the moment label showers by adding 1000 to ID

    fResHolder.evalRatio = evalRatio;
    fResHolder.concentration = concentration;
    fResHolder.coreHaloRatio = coreHaloRatio;
    fResHolder.conicalness = conicalness;
    fResHolder.dEdxStart = dEdxStart;
    fResHolder.dEdxEnd = dEdxEnd;
    if (dEdxPenultimate < 0.1)
      fResHolder.dEdxEndRatio = 1.0;
    else
      fResHolder.dEdxEndRatio = dEdxEnd / dEdxPenultimate;
    fResHolder.length = sortedObj.length;

    for (auto methodIter = fMVAMethods.begin(); methodIter != fMVAMethods.end(); ++methodIter) {
      fResHolder.mvaOutput[*methodIter] = fReader.EvaluateMVA(*methodIter);
    }
    result.push_back(fResHolder);
    //    util::CreateAssn(evt, result, *showerIter, showerAssns);
  }
}

//**********************************************************
void mvapid::MVAAlg::PrepareEvent(const AnaEventB& evt, const detinfo::DetectorClocksData& clockData){
//**********************************************************
  
  fHits.clear();
  fSpacePoints.clear();
  fTracks.clear();
  fShowers.clear();
  //  fSpacePointsToHits.clear();
  //  fHitsToSpacePoints.clear();
  fTracksToHits.clear();
  fTracksToSpacePoints.clear();
  fShowersToHits.clear();
  fShowersToSpacePoints.clear();

  fEventT0 = calo::trigger_offset(clockData);

  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<const AnaEventB*>(&evt)->Particles;
  int nParts           = static_cast<const AnaEventB*>(&evt)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);

    if (part->Hits[2].size() ==0) continue;
    
    Track* track = new Track(*part);
    fTracks.push_back(track);

    for (UInt_t j = 0; j < part->Hits[2].size(); j++){
      SpacePoint* spacePoint = new SpacePoint(part->HitPosition[2][j]);      
      fSpacePoints.push_back(spacePoint);
      fTracksToSpacePoints[track].push_back(spacePoint);
      AnaHitPD* hit = new AnaHitPD(part->Hits[2][j]);
      fHits.push_back(hit);
      fHitsToSpacePoints[hit] = spacePoint;
      fTracksToHits[track].push_back(hit);
    }
    /*
    Shower* shower = new Shower(*part);
    fShowers.push_back(shower);

    for (UInt_t j = 0; j < part->Hits[2].size(); j++){
      SpacePoint* spacePoint = new SpacePoint(part->HitPosition[2][j]);      
      fShowersToSpacePoints[shower].push_back(spacePoint);
      AnaHitPD& hit = part->Hits[2][j];
      fHitsToSpacePoints[&hit] = spacePoint;
      fShowersToHits[shower].push_back(&hit);
    }
    */
  }

  fVertex4Vect = anaUtils::ArrayToTLorentzVector(evt.TrueParticles[0]->Position);

}


//**********************************************************
void mvapid::MVAAlg::FitAndSortTrack(Track* track,
                                     int& isStoppingReco,
                                     mvapid::MVAAlg::SortedObj& sortedTrack){
//**********************************************************
  
  sortedTrack.hitMap.clear();
  TVector3 trackPoint, trackDir;
  this->LinFit(track, trackPoint, trackDir);

  TVector3 nearestPointStart, nearestPointEnd;

  //For single-particle events can opt to cheat vertex from start of primary trajectory.
  //Ok since in real events it should be possible to identify the true vertex.
  if (fCheatVertex) {
    if ((track->End<TVector3>() - fVertex4Vect.Vect()).Mag() >
        (track->Vertex<TVector3>() - fVertex4Vect.Vect()).Mag()) {
      nearestPointStart =
        trackPoint +
        trackDir * (trackDir.Dot(track->Vertex<TVector3>() - trackPoint) / trackDir.Mag2());
      nearestPointEnd = trackPoint + trackDir * (trackDir.Dot(track->End<TVector3>() - trackPoint) /
                                                 trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->End<TVector3>());
    }
    else {
      nearestPointStart =
        trackPoint +
        trackDir * (trackDir.Dot(track->End<TVector3>() - trackPoint) / trackDir.Mag2());
      nearestPointEnd =
        trackPoint +
        trackDir * (trackDir.Dot(track->Vertex<TVector3>() - trackPoint) / trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->Vertex<TVector3>());
      trackDir *= -1.;
    }
  }
  else {
    if (track->End<TVector3>().Z() >=
        track->Vertex<TVector3>().Z()) { //Otherwise assume particle is forward-going for now...
      nearestPointStart =
        trackPoint +
        trackDir * (trackDir.Dot(track->Vertex<TVector3>() - trackPoint) / trackDir.Mag2());
      nearestPointEnd = trackPoint + trackDir * (trackDir.Dot(track->End<TVector3>() - trackPoint) /
                                                 trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->End<TVector3>());
    }
    else {
      nearestPointStart =
        trackPoint +
        trackDir * (trackDir.Dot(track->End<TVector3>() - trackPoint) / trackDir.Mag2());
      nearestPointEnd =
        trackPoint +
        trackDir * (trackDir.Dot(track->Vertex<TVector3>() - trackPoint) / trackDir.Mag2());
      isStoppingReco = this->IsInActiveVol(track->Vertex<TVector3>());
    }

    if (trackDir.Z() <= 0) {
      trackDir.SetX(-trackDir.X());
      trackDir.SetY(-trackDir.Y());
      trackDir.SetZ(-trackDir.Z());
    }
  }

  sortedTrack.start = nearestPointStart;
  sortedTrack.end = nearestPointEnd;
  sortedTrack.dir = trackDir;
  sortedTrack.length = (nearestPointEnd - nearestPointStart).Mag();

  std::vector<AnaHitPD*> hits = fTracksToHits[track];

  for (auto hitIter = hits.begin(); hitIter != hits.end(); ++hitIter) {

    if (!fHitsToSpacePoints.count(*hitIter)) continue;
    SpacePoint* sp = fHitsToSpacePoints.at(*hitIter);

    TVector3 nearestPoint =
      trackPoint + trackDir * (trackDir.Dot(TVector3(sp->XYZ()) - trackPoint) / trackDir.Mag2());
    double lengthAlongTrack = (nearestPointStart - nearestPoint).Mag();
    sortedTrack.hitMap.insert(std::pair<double, AnaHitPD*>(lengthAlongTrack, *hitIter));
  }
}

//**********************************************************
void mvapid::MVAAlg::SortShower(Shower* shower,
                                int& isStoppingReco,
                                mvapid::MVAAlg::SortedObj& sortedShower){
//**********************************************************

  sortedShower.hitMap.clear();

  std::vector<AnaHitPD*> hits = fShowersToHits[shower];

  TVector3 showerEnd(0, 0, 0);
  double furthestHitFromStart = -999.9;
  for (auto hitIter = hits.begin(); hitIter != hits.end(); ++hitIter) {

    if (!fHitsToSpacePoints.count(*hitIter)) continue;
    SpacePoint* sp = fHitsToSpacePoints.at(*hitIter);
    if ((TVector3(sp->XYZ()) - shower->ShowerStart()).Mag() > furthestHitFromStart) {
      showerEnd = TVector3(sp->XYZ());
      furthestHitFromStart = (TVector3(sp->XYZ()) - shower->ShowerStart()).Mag();
    }
  }

  TVector3 showerPoint, showerDir;
  this->LinFitShower(shower, showerPoint, showerDir);

  TVector3 nearestPointStart, nearestPointEnd;

  //Ensure that shower is fitted in correct direction (assuming for now that particle moves in +z direction)

  if (fCheatVertex) {
    if ((showerEnd - fVertex4Vect.Vect()).Mag() >
        (shower->ShowerStart() - fVertex4Vect.Vect()).Mag()) {
      nearestPointStart =
        showerPoint +
        showerDir * (showerDir.Dot(shower->ShowerStart() - showerPoint) / showerDir.Mag2());
      nearestPointEnd =
        showerPoint + showerDir * (showerDir.Dot(showerEnd - showerPoint) / showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(showerEnd);
    }
    else {
      nearestPointStart =
        showerPoint + showerDir * (showerDir.Dot(showerEnd - showerPoint) / showerDir.Mag2());
      nearestPointEnd =
        showerPoint +
        showerDir * (showerDir.Dot(shower->ShowerStart() - showerPoint) / showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(shower->ShowerStart());
      showerDir *= -1.;
    }
  }
  else {
    if (showerEnd.Z() >= shower->ShowerStart().Z()) {
      nearestPointStart =
        showerPoint +
        showerDir * (showerDir.Dot(shower->ShowerStart() - showerPoint) / showerDir.Mag2());
      nearestPointEnd =
        showerPoint + showerDir * (showerDir.Dot(showerEnd - showerPoint) / showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(showerEnd);
    }
    else {
      nearestPointStart =
        showerPoint + showerDir * (showerDir.Dot(showerEnd - showerPoint) / showerDir.Mag2());
      nearestPointEnd =
        showerPoint +
        showerDir * (showerDir.Dot(shower->ShowerStart() - showerPoint) / showerDir.Mag2());
      isStoppingReco = this->IsInActiveVol(shower->ShowerStart());
    }

    if (showerDir.Z() <= 0) {
      showerDir.SetX(-showerDir.X());
      showerDir.SetY(-showerDir.Y());
      showerDir.SetZ(-showerDir.Z());
    }
  }

  sortedShower.start = nearestPointStart;
  sortedShower.end = nearestPointEnd;
  sortedShower.dir = showerDir;
  sortedShower.length = (nearestPointEnd - nearestPointStart).Mag();

  for (auto hitIter = hits.begin(); hitIter != hits.end(); ++hitIter) {

    if (!fHitsToSpacePoints.count(*hitIter)) continue;
    SpacePoint* sp = fHitsToSpacePoints.at(*hitIter);

    TVector3 nearestPoint =
      showerPoint +
      showerDir * (showerDir.Dot(TVector3(sp->XYZ()) - showerPoint) / showerDir.Mag2());
    double lengthAlongShower = (nearestPointStart - nearestPoint).Mag();
    sortedShower.hitMap.insert(
      std::pair<double, AnaHitPD*>(lengthAlongShower, *hitIter));
  }
}

//**********************************************************
void mvapid::MVAAlg::RunPCA(std::vector<AnaHitPD*>& hits,
                            std::vector<double>& eVals,
                            std::vector<double>& eVecs){
//**********************************************************

  TPrincipal* principal = new TPrincipal(3, "D");

  for (auto hitIter = hits.begin(); hitIter != hits.end(); ++hitIter) {

    if (fHitsToSpacePoints.count(*hitIter)) {
      principal->AddRow(fHitsToSpacePoints.at(*hitIter)->XYZ());
    }
  }

  // PERFORM PCA
  principal->MakePrincipals();
  // GET EIGENVALUES AND EIGENVECTORS
  for (unsigned int i = 0; i < 3; ++i) {
    eVals.push_back(principal->GetEigenValues()->GetMatrixArray()[i]);
  }

  for (unsigned int i = 0; i < 9; ++i) {
    eVecs.push_back(principal->GetEigenVectors()->GetMatrixArray()[i]);
  }
}

//**********************************************************
void mvapid::MVAAlg::_Var_Shape(const mvapid::MVAAlg::SortedObj& track,
                                double& coreHaloRatio,
                                double& concentration,
                                double& conicalness){
//**********************************************************  

  static const unsigned int conMinHits = 3;
  static const double minCharge = 0.1;
  static const double conFracRange = 0.2;
  static const double MoliereRadius = 10.1;
  static const double MoliereRadiusFraction = 0.2;

  double totalCharge = 0;
  double totalChargeStart = 0;
  double totalChargeEnd = 0;

  double chargeCore = 0;
  double chargeHalo = 0;
  double chargeCon = 0;
  unsigned int nHits = 0;

  //stuff for conicalness
  double chargeConStart = 0;
  double chargeConEnd = 0;
  unsigned int nHitsConStart = 0;
  unsigned int nHitsConEnd = 0;

  for (auto hitIter = track.hitMap.begin(); hitIter != track.hitMap.end(); ++hitIter) {
    if (fHitsToSpacePoints.count(hitIter->second)) {
      SpacePoint* sp = fHitsToSpacePoints.at(hitIter->second);

      double distFromTrackFit = ((TVector3(sp->XYZ()) - track.start).Cross(track.dir)).Mag();

      ++nHits;

      if (distFromTrackFit < MoliereRadiusFraction * MoliereRadius)
        chargeCore += hitIter->second->Integral();
      else
        chargeHalo += hitIter->second->Integral();

      totalCharge += hitIter->second->Integral();

      chargeCon += hitIter->second->Integral() / std::max(1.E-2, distFromTrackFit);
      if (hitIter->first / track.length < conFracRange) {
        chargeConStart += distFromTrackFit * distFromTrackFit * hitIter->second->Integral();
        ++nHitsConStart;
        totalChargeStart += hitIter->second->Integral();
      }
      else if (1. - hitIter->first / track.length < conFracRange) {
        chargeConEnd += distFromTrackFit * distFromTrackFit * hitIter->second->Integral();
        ++nHitsConEnd;
        totalChargeEnd += hitIter->second->Integral();
      }
    }
  }

  coreHaloRatio = chargeHalo / TMath::Max(1.0E-3, chargeCore);
  coreHaloRatio = TMath::Min(100.0, coreHaloRatio);
  concentration = chargeCon / totalCharge;
  if (nHitsConStart >= conMinHits && nHitsConEnd >= conMinHits && totalChargeEnd > minCharge &&
      sqrt(chargeConStart) > minCharge && totalChargeStart > minCharge) {
    conicalness = (sqrt(chargeConEnd) / totalChargeEnd) / (sqrt(chargeConStart) / totalChargeStart);
  }
  else {
    conicalness = 1.;
  }
}

//**********************************************************
double mvapid::MVAAlg::CalcSegmentdEdxFrac(const detinfo::DetectorClocksData& clock_data,
                                           const detinfo::DetectorPropertiesData& det_prop,
                                           const mvapid::MVAAlg::SortedObj& track,
                                           double start,
                                           double end){
//**********************************************************
  

  double trackLength = (track.end - track.start).Mag();
  return CalcSegmentdEdxDist(clock_data, det_prop, track, start * trackLength, end * trackLength);
}

//**********************************************************
double mvapid::MVAAlg::CalcSegmentdEdxDistAtEnd(const detinfo::DetectorClocksData& clock_data,
                                                const detinfo::DetectorPropertiesData& det_prop,
                                                const mvapid::MVAAlg::SortedObj& track,
                                                double distAtEnd){
//**********************************************************
  
  double trackLength = (track.end - track.start).Mag();
  return CalcSegmentdEdxDist(clock_data, det_prop, track, trackLength - distAtEnd, trackLength);
}

//**********************************************************
double mvapid::MVAAlg::CalcSegmentdEdxDist(const detinfo::DetectorClocksData& clock_data,
                                           const detinfo::DetectorPropertiesData& det_prop,
                                           const mvapid::MVAAlg::SortedObj& track,
                                           double start,
                                           double end){
//**********************************************************  


//  art::ServiceHandle<geo::Geometry const> geom;

  double totaldEdx = 0;
  unsigned int nHits = 0;

  // TODO: Harcoded wire pitch
  double wirePitch = 0.4792; 

  //Loop over hits again to calculate average dE/dx and shape variables
  for (auto hitIter = track.hitMap.begin(); hitIter != track.hitMap.end(); ++hitIter) {

    if (hitIter->first < start) continue;
    if (hitIter->first >= end) break;

    AnaHitPD* hit = hitIter->second;


    //Pitch to use in dEdx calculation
    double yzPitch = wirePitch;   // TODO
      //      geom->WirePitch(hit->WireID().Plane,
      //                      hit->WireID().TPC); //pitch not taking into account angle of track or shower
    double xComponent, pitch3D;

    TVector3 dir = track.dir;

    //This assumes equal numbers of TPCs in each cryostat and equal numbers of planes in each TPC
    int planeKey = 0;//hit->WireID().Cryostat * geom->NTPC(0) * geom->Nplanes(0, 0) +
                   //hit->WireID().TPC * geom->Nplanes(0, 0) + hit->WireID().Plane;

    if (fNormToWiresY.count(planeKey) && fNormToWiresZ.count(planeKey)) {
      TVector3 normToWires(0.0, fNormToWiresY.at(planeKey), fNormToWiresZ.at(planeKey));
      yzPitch = wirePitch/fabs(dir.Dot(normToWires));
        //        geom->WirePitch(hit->WireID().Plane, hit->WireID().TPC) / fabs(dir.Dot(normToWires));
    }

    xComponent = yzPitch * dir[0] / sqrt(dir[1] * dir[1] + dir[2] * dir[2]);
    pitch3D = sqrt(xComponent * xComponent + yzPitch * yzPitch);

    double dEdx = fCaloAlg.dEdx_AREA(clock_data, det_prop, *hit, pitch3D, fEventT0);
    if (dEdx < 50.) {
      ++nHits;
      totaldEdx += dEdx;
    }
  }

  return nHits ? totaldEdx / nHits : 0;

}

//**********************************************************
int mvapid::MVAAlg::LinFit( Track* track, TVector3& trackPoint, TVector3& trackDir){
//**********************************************************
  
  const std::vector<SpacePoint*>& sp = fTracksToSpacePoints.at(track);

  TGraph2D grFit(1);
  unsigned int iPt = 0;
  for (auto spIter = sp.begin(); spIter != sp.end(); ++spIter) {
    TVector3 point = (*spIter)->XYZ();
    grFit.SetPoint(iPt++, point.X(), point.Y(), point.Z());
  }

  //Lift from the ROOT line3Dfit.C tutorial
  ROOT::Fit::Fitter fitter;
  // make the functor object
  mvapid::MVAAlg::SumDistance2 sdist(&grFit);

  ROOT::Math::Functor fcn(sdist, 6);

  //Initial fit parameters from track start and end...
  TVector3 trackStart = track->Vertex<TVector3>();
  TVector3 trackEnd = track->End<TVector3>();
  trackDir = (trackEnd - trackStart).Unit();

  TVector3 x0 = trackStart - trackDir;
  TVector3 u = trackDir;

  double pStart[6] = {x0.X(), u.X(), x0.Y(), u.Y(), x0.Z(), u.Z()};

  fitter.SetFCN(fcn, pStart);

  bool ok = fitter.FitFCN();
  if (!ok) {
    trackPoint.SetXYZ(x0.X(), x0.Y(), x0.Z());
    trackDir.SetXYZ(u.X(), u.Y(), u.Z());
    trackDir = trackDir.Unit();
    return 1;
  }
  else {
    const ROOT::Fit::FitResult& result = fitter.Result();
    const double* parFit = result.GetParams();
    trackPoint.SetXYZ(parFit[0], parFit[2], parFit[4]);
    trackDir.SetXYZ(parFit[1], parFit[3], parFit[5]);
    trackDir = trackDir.Unit();
    return 0;
  }
}

//**********************************************************
int mvapid::MVAAlg::LinFitShower(Shower* shower,
                                 TVector3& showerPoint,
                                 TVector3& showerDir){
//**********************************************************  

  const std::vector<SpacePoint*>& sp = fShowersToSpacePoints.at(shower);

  TGraph2D grFit(1);
  unsigned int iPt = 0;
  for (auto spIter = sp.begin(); spIter != sp.end(); ++spIter) {
    TVector3 point = (*spIter)->XYZ();
    grFit.SetPoint(iPt++, point.X(), point.Y(), point.Z());
  }

  //Lift from the ROOT line3Dfit.C tutorial
  ROOT::Fit::Fitter fitter;
  // make the functor object
  mvapid::MVAAlg::SumDistance2 sdist(&grFit);

  ROOT::Math::Functor fcn(sdist, 6);

  //Initial fit parameters from shower start and end...
  TVector3 showerStart = shower->ShowerStart();
  showerDir = shower->Direction().Unit();

  TVector3 x0 = showerStart - showerDir;
  TVector3 u = showerDir;

  double pStart[6] = {x0.X(), u.X(), x0.Y(), u.Y(), x0.Z(), u.Z()};

  fitter.SetFCN(fcn, pStart);

  bool ok = fitter.FitFCN();
  if (!ok) {
    showerPoint.SetXYZ(x0.X(), x0.Y(), x0.Z());
    showerDir.SetXYZ(u.X(), u.Y(), u.Z());
    showerDir = showerDir.Unit();
    return 1;
  }
  else {
    const ROOT::Fit::FitResult& result = fitter.Result();
    const double* parFit = result.GetParams();
    showerPoint.SetXYZ(parFit[0], parFit[2], parFit[4]);
    showerDir.SetXYZ(parFit[1], parFit[3], parFit[5]);
    showerDir = showerDir.Unit();
    return 0;
  }
}
