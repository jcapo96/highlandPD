#include "pdNeutralUtils.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdCreationUtils.hxx"
#include "pdNeutralHelpers.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>
#include <tuple>

namespace pdNeutralUtils {

//***************************************************************
std::pair<Int_t, Int_t> pdNeutralUtils::CalculateNeutralScore(
    AnaNeutralParticlePD* neutralParticle,
    AnaVertexPD* vertex,
    AnaParticlePD* parentParticle,
    AnaEventB& event,
    const std::unordered_map<Int_t, AnaParticlePD*>& particleByUniqueID) {
//***************************************************************

  // Initialize return values
  Int_t NPotentialParents = 0;
  Int_t NRecoHitsInVertex = 0;

  if (!neutralParticle || !vertex || !parentParticle) {
    return std::make_pair(NPotentialParents, NRecoHitsInVertex);
  }

  // Get array of particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // Get parameters
  double daughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
  double cylinderRadius = ND::params().GetParameterD("neutralKaonAnalysis.CylinderRadius");

  TVector3 vertexPos(vertex->Position[0], vertex->Position[1], vertex->Position[2]);

  // === ENHANCED NEUTRAL PARTICLE SCORING ALGORITHM ===
  // LOWER score = more neutral-like (few hits, far from path, not aligned)
  // HIGHER score = more charged-like (many hits, close to path, aligned)

  // Define neutral particle trajectory for scoring. Use Pandora position when valid so that
  // the same candidates are kept by FilterNeutralsByScore regardless of UsePandora (resolution studies).
  TVector3 neutralStart(neutralParticle->PositionStart[0],
                        neutralParticle->PositionStart[1],
                        neutralParticle->PositionStart[2]);
  TVector3 neutralEnd(neutralParticle->PositionEnd[0],
                      neutralParticle->PositionEnd[1],
                      neutralParticle->PositionEnd[2]);
  if (vertex->PositionPandora[0] > -900.f && vertex->PositionPandora[1] > -900.f && vertex->PositionPandora[2] > -900.f) {
    neutralEnd.SetXYZ(vertex->PositionPandora[0], vertex->PositionPandora[1], vertex->PositionPandora[2]);
  }
  TVector3 neutralDirection = (neutralEnd - neutralStart).Unit();
  double neutralLength = (neutralEnd - neutralStart).Mag();

  // Collect all hits within the cylinder for analysis
  std::vector<TVector3> hitsInCylinder;
  std::vector<double> perpDistances;
  std::vector<double> longitudinalProjections;
  double totalDistance = 0.0;
  int nHitsInCylinder = 0;

  for(int j = 0; j < nParts; j++){
    AnaParticlePD* potentialParent = static_cast<AnaParticlePD*>(parts[j]);
    if(!potentialParent) continue;

    // Skip vertex particles and current parent
    if(vertex->Particles.size() >= 2){
      if(potentialParent->UniqueID == vertex->Particles[0]->UniqueID ||
         potentialParent->UniqueID == vertex->Particles[1]->UniqueID ||
         potentialParent->UniqueID == parentParticle->UniqueID) continue;
    }

    TVector3 potentialParentEnd(potentialParent->PositionEnd[0],
                                 potentialParent->PositionEnd[1],
                                 potentialParent->PositionEnd[2]);
    double potentialParentDistance = (potentialParentEnd - vertexPos).Mag();

    if(potentialParentDistance < daughterDistance){
      // OPTIMIZATION: Count potential parents (merged from CreateNeutral loop)
      NPotentialParents++;

      // Check hits from this particle
      for(int plane = 0; plane < 3; plane++){
        for(size_t k = 0; k < potentialParent->Hits[plane].size(); k++){
          AnaHitPD& hit = potentialParent->Hits[plane][k];
          if(hit.Position.Z() == -999) continue;

          TVector3 hitPos = hit.Position;

          // Check if hit is within cylinder
          TVector3 startToHit = hitPos - neutralStart;
          double projection = startToHit.Dot(neutralDirection);

          // Only consider hits along the neutral particle path
          if(projection >= 0 && projection <= neutralLength){
            TVector3 projectionPoint = neutralStart + projection * neutralDirection;
            double perpDistance = (hitPos - projectionPoint).Mag();

            if(perpDistance < cylinderRadius){
              nHitsInCylinder++;
              hitsInCylinder.push_back(hitPos);
              perpDistances.push_back(perpDistance);
              longitudinalProjections.push_back(projection);
              totalDistance += perpDistance;
              // OPTIMIZATION: Count hits in vertex (merged from CreateNeutral loop)
              NRecoHitsInVertex++;
            }
          }
        }
      }
    }
  }

  // === 1. ENHANCED DISTANCE METRICS ===

  // Calculate average perpendicular distance
  double avgDistance = (nHitsInCylinder > 0) ? totalDistance / nHitsInCylinder : 0.0;

  // Calculate RMS of perpendicular distances
  double rmsDistance = 0.0;
  if(nHitsInCylinder > 1){
    double sumSquaredDiff = 0.0;
    for(double dist : perpDistances){
      sumSquaredDiff += (dist - avgDistance) * (dist - avgDistance);
    }
    rmsDistance = TMath::Sqrt(sumSquaredDiff / nHitsInCylinder);
  }

  // === 2. LONGITUDINAL SPAN CALCULATION ===

  double longitudinalSpan = 0.0;
  if(nHitsInCylinder > 0 && neutralLength > 0){
    double minProj = *std::min_element(longitudinalProjections.begin(), longitudinalProjections.end());
    double maxProj = *std::max_element(longitudinalProjections.begin(), longitudinalProjections.end());
    double span = maxProj - minProj;
    longitudinalSpan = span / neutralLength;  // Fraction of path covered by hits
  }

  // === 3. ROBUST ALIGNMENT CALCULATION ===

  double hitsAlignment = 0.0;  // Default: no alignment

  if(nHitsInCylinder == 1){
    // Single hit: check if vector from start to hit is aligned with neutral direction
    TVector3 hitVector = (hitsInCylinder[0] - neutralStart);
    if(hitVector.Mag() > 0){
      hitVector = hitVector.Unit();
      hitsAlignment = TMath::Abs(hitVector.Dot(neutralDirection));
    }
  }
  else if(nHitsInCylinder == 2){
    // Two hits: fit line between them
    TVector3 hitLine = (hitsInCylinder[1] - hitsInCylinder[0]);
    if(hitLine.Mag() > 0){
      hitLine = hitLine.Unit();
      hitsAlignment = TMath::Abs(hitLine.Dot(neutralDirection));
    }
  }
  else if(nHitsInCylinder >= 3){
    // Three or more hits: use PCA
    // Calculate centroid of hits
    TVector3 centroid(0, 0, 0);
    for(const auto& hit : hitsInCylinder){
      centroid += hit;
    }
    centroid *= (1.0 / hitsInCylinder.size());

    // Build covariance matrix
    double cov_xx = 0, cov_yy = 0, cov_zz = 0;
    double cov_xy = 0, cov_xz = 0, cov_yz = 0;

    for(const auto& hit : hitsInCylinder){
      TVector3 diff = hit - centroid;
      cov_xx += diff.X() * diff.X();
      cov_yy += diff.Y() * diff.Y();
      cov_zz += diff.Z() * diff.Z();
      cov_xy += diff.X() * diff.Y();
      cov_xz += diff.X() * diff.Z();
      cov_yz += diff.Y() * diff.Z();
    }

    // Normalize
    double norm = 1.0 / hitsInCylinder.size();
    cov_xx *= norm;
    cov_yy *= norm;
    cov_zz *= norm;
    cov_xy *= norm;
    cov_xz *= norm;
    cov_yz *= norm;

    // Find principal eigenvector (direction of maximum variance)
    // Use power iteration method
    TVector3 fittedDirection(1, 1, 1);  // Initial guess
    for(int iter = 0; iter < 20; iter++){
      TVector3 newDir;
      newDir.SetX(cov_xx * fittedDirection.X() + cov_xy * fittedDirection.Y() + cov_xz * fittedDirection.Z());
      newDir.SetY(cov_xy * fittedDirection.X() + cov_yy * fittedDirection.Y() + cov_yz * fittedDirection.Z());
      newDir.SetZ(cov_xz * fittedDirection.X() + cov_yz * fittedDirection.Y() + cov_zz * fittedDirection.Z());

      double mag = newDir.Mag();
      if(mag > 0){
        fittedDirection = newDir * (1.0 / mag);
      }
    }

    // Calculate alignment as absolute value of dot product
    hitsAlignment = TMath::Abs(fittedDirection.Dot(neutralDirection));
  }

  // === 4. CALCULATE SCORE COMPONENTS ===

  // A) Hit density (hits per unit distance)
  // This will be a multiplicative factor for hit-dependent scores
  // Ensures particles with 0 hits get 0 score for hit-dependent component
  double hitDensity = 0.0;
  if(neutralLength > 0){
    hitDensity = nHitsInCylinder / neutralLength;
  }

  // B) Enhanced distance score - combines average and RMS
  // Scattered hits (high RMS) are more neutral-like
  // Concentrated hits (low RMS) are more charged-like
  double distanceScore = 0.0;
  if(nHitsInCylinder > 0 && cylinderRadius > 0){
    double avgFactor = 1.0 - TMath::Min(avgDistance / cylinderRadius, 1.0);
    double rmsFactor = (nHitsInCylinder > 1) ? (1.0 - TMath::Min(rmsDistance / cylinderRadius, 1.0)) : avgFactor;
    // Lower RMS (concentrated) → higher rmsFactor → higher score (charged-like)
    distanceScore = avgFactor * rmsFactor * 100.0;
  }

  // C) Alignment score - already calculated (0 to 1, higher = more aligned)
  double alignmentScore = hitsAlignment * 100.0;

  // D) Clustering score - longitudinal span
  // High span (hits spread along path) → charged-like
  // Low span (localized hits) → could be neutral or random
  double clusterScore = longitudinalSpan * 100.0;

  // === 5. CALCULATE SIBLING ALIGNMENT SCORE ===
  // Look at alignment between all parent's daughters (including vertex particles) and neutral direction

  double siblingAlignmentScore = 0.0;
  int nSiblings = 0;
  double totalSiblingAlignment = 0.0;

  if(parentParticle && parentParticle->Daughters.size() > 0){
    // Loop through parent's daughters (including vertex particles)
    for(size_t i = 0; i < parentParticle->Daughters.size(); i++){
      AnaParticlePD* sibling = static_cast<AnaParticlePD*>(parentParticle->Daughters[i]);
      if(!sibling) continue;

      // Check if sibling has valid direction
      if(sibling->DirectionStart[0] < -900) continue;

      // Calculate alignment between sibling direction and neutral direction
      TVector3 siblingDir(sibling->DirectionStart[0],
                         sibling->DirectionStart[1],
                         sibling->DirectionStart[2]);
      if(siblingDir.Mag() > 0){
        siblingDir = siblingDir.Unit();
        double alignment = TMath::Abs(siblingDir.Dot(neutralDirection));
        totalSiblingAlignment += alignment;
        nSiblings++;
      }
    }

    // Average sibling alignment
    if(nSiblings > 0){
      siblingAlignmentScore = (totalSiblingAlignment / nSiblings) * 100.0;
    }
  }

  // === 6. COMBINED SCORING ===
  // neutralScore = (nhits/length) * (alignmentScore + distanceScore + clusterScore) + siblingAlignmentScore
  // Equal weights for hit-dependent components, multiplied by hit density
  // Sibling alignment is independent of hit density

  double hitDependentScore = hitDensity * (alignmentScore + distanceScore + clusterScore);
  double neutralScore = hitDependentScore + siblingAlignmentScore;

  // Store all metrics in neutral particle
  neutralParticle->NeutralScore = neutralScore;
  neutralParticle->HitsAlignment = hitsAlignment;
  neutralParticle->NHitsInCylinder = nHitsInCylinder;
  neutralParticle->HitsAvgDistance = avgDistance;
  neutralParticle->HitsRMSDistance = rmsDistance;
  neutralParticle->HitsLongitudinalSpan = longitudinalSpan;

  // OPTIMIZATION: Return counts that were previously calculated in a separate loop
  return std::make_pair(NPotentialParents, NRecoHitsInVertex);
}

//***************************************************************
AnaTrueEquivalentNeutralParticlePD* pdNeutralUtils::FillTrueEquivalentNeutralParticle(
    AnaVertexPD* vertex,
    AnaParticlePD* parentParticle) {
//***************************************************************

  if (!vertex || !parentParticle) {
    return nullptr;
  }

  // Create true equivalent neutral particle
  AnaTrueEquivalentNeutralParticlePD* trueEquivalentNeutralParticle = new AnaTrueEquivalentNeutralParticlePD();
  trueEquivalentNeutralParticle->TrueEquivalentVertex = static_cast<AnaTrueEquivalentVertexPD*>(vertex->TrueEquivalentVertex);
  trueEquivalentNeutralParticle->TrueParent = static_cast<AnaTrueParticlePD*>(parentParticle->TrueObject);

  // Add the position and direction of the true particle
  if (parentParticle->TrueObject) {
    AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(parentParticle->TrueObject);
    trueEquivalentNeutralParticle->Position[0] = trueParticle->PositionEnd[0];
    trueEquivalentNeutralParticle->Position[1] = trueParticle->PositionEnd[1];
    trueEquivalentNeutralParticle->Position[2] = trueParticle->PositionEnd[2];
  } else {
    trueEquivalentNeutralParticle->Position[0] = -999.0;
    trueEquivalentNeutralParticle->Position[1] = -999.0;
    trueEquivalentNeutralParticle->Position[2] = -999.0;
  }

  // Add the end position of the true vertex
  if (vertex->TrueEquivalentVertex) {
    trueEquivalentNeutralParticle->PositionEnd[0] = vertex->TrueEquivalentVertex->Position[0];
    trueEquivalentNeutralParticle->PositionEnd[1] = vertex->TrueEquivalentVertex->Position[1];
    trueEquivalentNeutralParticle->PositionEnd[2] = vertex->TrueEquivalentVertex->Position[2];
  } else {
    trueEquivalentNeutralParticle->PositionEnd[0] = -999.0;
    trueEquivalentNeutralParticle->PositionEnd[1] = -999.0;
    trueEquivalentNeutralParticle->PositionEnd[2] = -999.0;
  }

  // Calculate the direction of the true particle
  Float_t trueDirectionStart[3];
  trueDirectionStart[0] = trueEquivalentNeutralParticle->PositionEnd[0] - trueEquivalentNeutralParticle->Position[0];
  trueDirectionStart[1] = trueEquivalentNeutralParticle->PositionEnd[1] - trueEquivalentNeutralParticle->Position[1];
  trueDirectionStart[2] = trueEquivalentNeutralParticle->PositionEnd[2] - trueEquivalentNeutralParticle->Position[2];

  Float_t trueNorm = sqrt(trueDirectionStart[0]*trueDirectionStart[0] +
                          trueDirectionStart[1]*trueDirectionStart[1] +
                          trueDirectionStart[2]*trueDirectionStart[2]);
  if (trueNorm > 0) {
    trueDirectionStart[0] /= trueNorm;
    trueDirectionStart[1] /= trueNorm;
    trueDirectionStart[2] /= trueNorm;
  } else {
    trueDirectionStart[0] = -999.0;
    trueDirectionStart[1] = -999.0;
    trueDirectionStart[2] = -999.0;
  }

  trueEquivalentNeutralParticle->Direction[0] = trueDirectionStart[0];
  trueEquivalentNeutralParticle->Direction[1] = trueDirectionStart[1];
  trueEquivalentNeutralParticle->Direction[2] = trueDirectionStart[2];

  Float_t trueDirectionEnd[3];
  trueDirectionEnd[0] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[0] : -999;
  trueDirectionEnd[1] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[1] : -999;
  trueDirectionEnd[2] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[2] : -999;

  Float_t trueNormEnd = sqrt(trueDirectionEnd[0]*trueDirectionEnd[0] +
                             trueDirectionEnd[1]*trueDirectionEnd[1] +
                             trueDirectionEnd[2]*trueDirectionEnd[2]);
  if (trueNormEnd > 0) {
    trueDirectionEnd[0] /= trueNormEnd;
    trueDirectionEnd[1] /= trueNormEnd;
    trueDirectionEnd[2] /= trueNormEnd;
  } else {
    trueDirectionEnd[0] = -999.0;
    trueDirectionEnd[1] = -999.0;
    trueDirectionEnd[2] = -999.0;
  }

  trueEquivalentNeutralParticle->DirectionEnd[0] = trueDirectionEnd[0];
  trueEquivalentNeutralParticle->DirectionEnd[1] = trueDirectionEnd[1];
  trueEquivalentNeutralParticle->DirectionEnd[2] = trueDirectionEnd[2];

  Float_t trueLength = sqrt(pow(trueEquivalentNeutralParticle->PositionEnd[0]-trueEquivalentNeutralParticle->Position[0],2)+
                            pow(trueEquivalentNeutralParticle->PositionEnd[1]-trueEquivalentNeutralParticle->Position[1],2)+
                            pow(trueEquivalentNeutralParticle->PositionEnd[2]-trueEquivalentNeutralParticle->Position[2],2));
  trueEquivalentNeutralParticle->Length = trueLength;

  // Calculate invariant mass and momentum for true equivalent neutral particle
  Float_t trueInvariantMass = -999;
  Float_t trueMomentumEnd = -999;

  if (vertex->TrueEquivalentVertex && vertex->TrueEquivalentVertex->TrueParticles.size() >= 2) {
    AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[0]);
    AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[1]);

    if (trueParticle1 && trueParticle2 &&
        trueParticle1->Momentum > 0 && trueParticle2->Momentum > 0 &&
        trueParticle1->Momentum != -999 && trueParticle2->Momentum != -999) {
      const Float_t pionMass = 0.13957;
      trueInvariantMass = pdAnaUtils::ComputeTrueInvariantMass(*trueParticle1, *trueParticle2, pionMass, pionMass);

      // Calculate total momentum (sum of daughter momenta)
      trueMomentumEnd = trueParticle1->Momentum + trueParticle2->Momentum;
    }
  }

  trueEquivalentNeutralParticle->Mass = trueInvariantMass;
  trueEquivalentNeutralParticle->MomentumEnd = trueMomentumEnd;

  return trueEquivalentNeutralParticle;
}

//***************************************************************
void pdNeutralUtils::CalculateVertexDegeneracies(AnaEventB& event,
                                                   AnaCreationVertexPD* creationVertex,
                                                   AnaAnnihilationVertexPD* annihilationVertex){
//***************************************************************

  // Calculate creation vertex degeneracy
  // Count particles whose start position is within the creation vertex radius, excluding annihilation vertex particles
  int creationDegeneracy = 0;
  std::vector<double> creationDistances; // Store distances for DegDist array
  double creationVertexRadius = ND::params().GetParameterD("neutralKaonAnalysis.CreationVertexRadius");
  TVector3 creationPos(creationVertex->Position[0], creationVertex->Position[1], creationVertex->Position[2]);

  // Initialize DegDist array
  for (int i = 0; i < 5; i++) {
    creationVertex->DegDist[i] = -999.0;
  }

  for(int p = 0; p < event.nParticles; p++){
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[p]);
    if(!particle) continue;

    // Skip particles in the annihilation vertex
    bool isInAnnihilationVertex = false;
    for (int i = 0; i < annihilationVertex->NParticles; i++) {
      if (annihilationVertex->Particles[i] && particle->UniqueID == annihilationVertex->Particles[i]->UniqueID) {
        isInAnnihilationVertex = true;
        break;
      }
    }
    if (isInAnnihilationVertex) continue;

    // Calculate distance from particle start to creation vertex
    TVector3 particleStart(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
    double distance = (particleStart - creationPos).Mag();

    // Count if within radius and store distance
    if(distance < creationVertexRadius) {
      creationDegeneracy++;
      creationDistances.push_back(distance);
    }
  }

  creationVertex->Degeneracy = creationDegeneracy;

  // Sort distances and fill DegDist array (up to 5)
  std::sort(creationDistances.begin(), creationDistances.end());
  for (size_t i = 0; i < creationDistances.size() && i < 5; i++) {
    creationVertex->DegDist[i] = static_cast<Float_t>(creationDistances[i]);
  }

  // Calculate annihilation vertex degeneracy and distances
  // Count particles whose start position is within the annihilation vertex radius, excluding creation vertex particles
  int annihilationDegeneracy = 0;
  std::vector<double> annihilationDistances; // Store distances for DegDist array
  double maxDaughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");
  TVector3 annihilationPos(annihilationVertex->Position[0], annihilationVertex->Position[1], annihilationVertex->Position[2]);

  // Initialize DegDist array
  for (int i = 0; i < 5; i++) {
    annihilationVertex->DegDist[i] = -999.0;
  }

  for(int p = 0; p < event.nParticles; p++){
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(event.Particles[p]);
    if(!particle) continue;

    // Skip the creation vertex secondary particle (if it exists)
    if(creationVertex->SecondParticle && particle->UniqueID == creationVertex->SecondParticle->UniqueID) continue;

    // Skip the beam particle
    if(particle->UniqueID == creationVertex->BeamParticle->UniqueID) continue;

    // Calculate distance from particle start to annihilation vertex
    TVector3 particleStart(particle->PositionStart[0], particle->PositionStart[1], particle->PositionStart[2]);
    double distance = (particleStart - annihilationPos).Mag();

    // Count if within radius and store distance
    if(distance < maxDaughterDistance) {
      annihilationDegeneracy++;
      annihilationDistances.push_back(distance);
    }
  }

  annihilationVertex->Degeneracy = annihilationDegeneracy;

  // Sort distances and fill DegDist array (up to 5)
  std::sort(annihilationDistances.begin(), annihilationDistances.end());
  for (size_t i = 0; i < annihilationDistances.size() && i < 5; i++) {
    annihilationVertex->DegDist[i] = static_cast<Float_t>(annihilationDistances[i]);
  }
}

//***************************************************************
std::tuple<double, double, double> pdNeutralUtils::CalculateRawScoreComponents(AnaNeutralParticlePD* neutralParticle){
//***************************************************************

  if (!neutralParticle || !neutralParticle->CreationVertex) {
    return std::make_tuple(999.0, 1.0, 999.0); // Worst possible scores
  }

  // 1. ProtonScore from CreationVertex (lower is better)
  double protonScore = 999.0;
  if (neutralParticle->CreationVertex) {
    protonScore = static_cast<double>(neutralParticle->CreationVertex->ProtonScore);
  }

  // 2. Alignment score: alignment between parent direction and neutral direction
  // Higher alignment = better, so we use (1 - alignment) as the score (lower is better)
  double alignmentScore = 1.0; // Worst case: no alignment
  if (neutralParticle->CreationVertex && neutralParticle->CreationVertex->BeamParticle) {
    AnaParticlePD* parentParticle = neutralParticle->CreationVertex->BeamParticle;

    // Get parent direction (use DirectionEnd for parent)
    TVector3 parentDir(parentParticle->DirectionEnd[0],
                       parentParticle->DirectionEnd[1],
                       parentParticle->DirectionEnd[2]);

    // Get neutral direction
    TVector3 neutralDir(neutralParticle->DirectionStart[0],
                        neutralParticle->DirectionStart[1],
                        neutralParticle->DirectionStart[2]);

    // Normalize directions
    if (parentDir.Mag() > 0 && neutralDir.Mag() > 0) {
      parentDir = parentDir.Unit();
      neutralDir = neutralDir.Unit();

      // Calculate alignment (dot product, ranges from -1 to 1)
      double alignment = TMath::Abs(parentDir.Dot(neutralDir));

      // Convert to score: (1 - alignment) so that higher alignment = lower score
      alignmentScore = 1.0 - alignment;
    }

  }

  // 3. Number of hits in cylinder (already calculated, lower is better)
  double hitsScore = static_cast<double>(neutralParticle->NHitsInCylinder);

  return std::make_tuple(protonScore, alignmentScore, hitsScore);
}

//***************************************************************
void pdNeutralUtils::NormalizeNeutralParticleScores(std::vector<AnaNeutralParticlePD*>& neutralParticles){
//***************************************************************

  if (neutralParticles.empty()) return;

  // Collect all alignment scores (dot product with beam direction)
  std::vector<double> alignmentScores;

  for (AnaNeutralParticlePD* neutralParticle : neutralParticles) {
    if (!neutralParticle) continue;

    // Calculate alignment: dot product between neutral direction and beam particle direction
    double alignment = -1.0; // Default to worst case
    if (neutralParticle->CreationVertex && neutralParticle->CreationVertex->BeamParticle) {
      AnaParticlePD* beamParticle = neutralParticle->CreationVertex->BeamParticle;

      // Get beam particle direction (use DirectionEnd for beam)
      TVector3 beamDir(beamParticle->DirectionEnd[0],
                       beamParticle->DirectionEnd[1],
                       beamParticle->DirectionEnd[2]);

      // Get neutral direction
      TVector3 neutralDir(neutralParticle->DirectionStart[0],
                          neutralParticle->DirectionStart[1],
                          neutralParticle->DirectionStart[2]);

      // Normalize directions
      if (beamDir.Mag() > 0 && neutralDir.Mag() > 0) {
        beamDir = beamDir.Unit();
        neutralDir = neutralDir.Unit();

        // Calculate alignment (dot product, ranges from -1 to 1)
        // Positive means forward along beam direction
        alignment = beamDir.Dot(neutralDir);
      }
    }
    alignmentScores.push_back(alignment);
  }

  if (alignmentScores.empty()) return;

  // Find min and max for alignment
  double alignmentMin = *std::min_element(alignmentScores.begin(), alignmentScores.end());
  double alignmentMax = *std::max_element(alignmentScores.begin(), alignmentScores.end());

  // Normalize alignment and calculate final score
  for (size_t i = 0; i < neutralParticles.size(); i++) {
    AnaNeutralParticlePD* neutralParticle = neutralParticles[i];
    if (!neutralParticle) continue;

    // Normalize AlignmentScore: best (highest, closest to 1) → 1.0, worst (lowest, closest to -1) → 0.0
    double normalizedAlignment = 0.5;
    if (alignmentMax > alignmentMin) {
      normalizedAlignment = (alignmentScores[i] - alignmentMin) / (alignmentMax - alignmentMin);
    }

    // Final score: use normalized alignment (higher is better, so we use 1 - normalized for lower-is-better sorting)
    // But we want higher alignment = lower score for sorting, so invert
    neutralParticle->NeutralScore = 1.0 - normalizedAlignment;
  }
}

//***************************************************************
std::vector<AnaNeutralParticlePD*> pdNeutralUtils::FilterNeutralsByScore(std::vector<AnaNeutralParticlePD*>& neutralParticles){
//***************************************************************

  // Filter neutral particles: keep only one per annihilation vertex (lowest NeutralScore)
  // Sort by NeutralScore (ascending - lower is better, more neutral-like)
  std::sort(neutralParticles.begin(), neutralParticles.end(),
    [](const AnaNeutralParticlePD* a, const AnaNeutralParticlePD* b) {
      // Handle invalid scores
      if (a->NeutralScore < -900 && b->NeutralScore < -900) return false;
      if (a->NeutralScore < -900) return false;  // Invalid score goes to back
      if (b->NeutralScore < -900) return true;   // Invalid score goes to back
      return a->NeutralScore < b->NeutralScore;
    });

  // Track which annihilation vertices have been used
  std::set<AnaAnnihilationVertexPD*> usedAnnihilationVertices;

  // Output vector for selected neutral particles
  std::vector<AnaNeutralParticlePD*> selectedNeutralParticles;
  std::set<AnaNeutralParticlePD*> selectedSet;  // For quick lookup

  // Iterate through sorted neutral particles (best scores first)
  // This ensures each annihilation vertex gets exactly one neutral particle (the best one)
  for (AnaNeutralParticlePD* neutralParticle : neutralParticles) {
    if (!neutralParticle || !neutralParticle->AnnihilationVertex) continue;

    AnaAnnihilationVertexPD* annihilationVtx = neutralParticle->AnnihilationVertex;

    // Check if this annihilation vertex is already used
    // Since particles are sorted by score (best first), the first one we encounter
    // for each annihilation vertex is the best one
    if (usedAnnihilationVertices.find(annihilationVtx) == usedAnnihilationVertices.end()) {
      // Annihilation vertex not used yet, so we keep this neutral particle (the best one for this vertex)
      selectedNeutralParticles.push_back(neutralParticle);
      selectedSet.insert(neutralParticle);

      // Mark annihilation vertex as used - no other neutral particle can use it
      usedAnnihilationVertices.insert(annihilationVtx);
    }
    // If annihilation vertex is already used, skip this neutral particle
    // (a better one was already selected)
  }

  // Delete neutral particles that were not selected
  for (AnaNeutralParticlePD* neutralParticle : neutralParticles) {
    if (selectedSet.find(neutralParticle) == selectedSet.end()) {
      delete neutralParticle;
    }
  }

  return selectedNeutralParticles;
}

//***************************************************************
std::vector<AnaNeutralParticlePD*> pdNeutralUtils::CreateNeutrals(AnaEventB& event,
                                                                  const std::vector<AnaCreationVertexPD*>& creationVertices,
                                                                  const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices){
//***************************************************************

  std::vector<AnaNeutralParticlePD*> allNeutralParticles;

  int neutralParticleID = 0;

  // OPTIMIZATION: Build hash map once for all neutral particles
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  std::unordered_map<Int_t, AnaParticlePD*> particleByUniqueID;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part){
      particleByUniqueID[part->UniqueID] = part;
    }
  }

  // Collect all particles from annihilation vertices to exclude from creation vertices
  std::set<Int_t> annihilationVertexParticleIDs;
  for (AnaAnnihilationVertexPD* annihilationVtx : annihilationVertices) {
    if (!annihilationVtx) continue;
    for (AnaParticlePD* particle : annihilationVtx->Particles) {
      if (particle) {
        annihilationVertexParticleIDs.insert(particle->UniqueID);
      }
    }
  }

  // Create neutral particles for ALL combinations of creation and annihilation vertices
  for(size_t a = 0; a < annihilationVertices.size(); a++){
    AnaAnnihilationVertexPD* annihilationVtx = annihilationVertices[a];
    if(!annihilationVtx) continue;

    for(size_t c = 0; c < creationVertices.size(); c++){
      AnaCreationVertexPD* creationVtx = creationVertices[c];
      if(!creationVtx) continue;

      // Skip if creation vertex's beam particle or second particle is already in an annihilation vertex
      if (creationVtx->BeamParticle &&
          annihilationVertexParticleIDs.find(creationVtx->BeamParticle->UniqueID) != annihilationVertexParticleIDs.end()) {
        continue;
      }
      if (creationVtx->SecondParticle &&
          annihilationVertexParticleIDs.find(creationVtx->SecondParticle->UniqueID) != annihilationVertexParticleIDs.end()) {
        continue;
      }

      // Sanity check: creation and annihilation vertices must be separated (K0 flight distance)
      // Reject pairs where they nearly coincide (would indicate track/vertex confusion)
      double minCreationAnnihilationSep = 5.0;
      if (ND::params().HasParameter("neutralKaonAnalysis.MinCreationAnnihilationSeparation")) {
        minCreationAnnihilationSep = ND::params().GetParameterD("neutralKaonAnalysis.MinCreationAnnihilationSeparation");
      }
      TVector3 creationPos(creationVtx->Position[0], creationVtx->Position[1], creationVtx->Position[2]);
      TVector3 annihilationPos(annihilationVtx->Position[0], annihilationVtx->Position[1], annihilationVtx->Position[2]);
      if ((creationPos - annihilationPos).Mag() < minCreationAnnihilationSep) {
        continue;
      }

      // Skip beam-only creation vertices (DistanceScore < 0.1): no proton constraint, poorly defined
      if (creationVtx->DistanceScore < 0.1) {
        continue;
      }

      // Create neutral particle for this combination
      AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();

      // Set vertices and parent
      neutralParticle->AnnihilationVertex = annihilationVtx;
      neutralParticle->CreationVertex = creationVtx;
      neutralParticle->Parent = creationVtx->BeamParticle;
      neutralParticle->UniqueID = neutralParticleID++;

      // Calculate degeneracies for both vertices
      CalculateVertexDegeneracies(event, creationVtx, annihilationVtx);

      // Set neutral particle start position from creation vertex
      neutralParticle->PositionStart[0] = creationVtx->Position[0];
      neutralParticle->PositionStart[1] = creationVtx->Position[1];
      neutralParticle->PositionStart[2] = creationVtx->Position[2];
      neutralParticle->PositionStart[3] = creationVtx->BeamParticle->PositionEnd[3];

      // Set end position from annihilation vertex
      neutralParticle->PositionEnd[0] = annihilationVtx->Position[0];
      neutralParticle->PositionEnd[1] = annihilationVtx->Position[1];
      neutralParticle->PositionEnd[2] = annihilationVtx->Position[2];
      neutralParticle->PositionEnd[3] = -999;

      // Calculate start direction
      Float_t directionStart[3];
      directionStart[0] = neutralParticle->PositionEnd[0] - neutralParticle->PositionStart[0];
      directionStart[1] = neutralParticle->PositionEnd[1] - neutralParticle->PositionStart[1];
      directionStart[2] = neutralParticle->PositionEnd[2] - neutralParticle->PositionStart[2];

      Float_t norm = sqrt(directionStart[0]*directionStart[0] +
                          directionStart[1]*directionStart[1] +
                          directionStart[2]*directionStart[2]);
      if (norm > 0) {
        directionStart[0] /= norm;
        directionStart[1] /= norm;
        directionStart[2] /= norm;
      }

      neutralParticle->DirectionStart[0] = directionStart[0];
      neutralParticle->DirectionStart[1] = directionStart[1];
      neutralParticle->DirectionStart[2] = directionStart[2];

      // Set end direction from annihilation vertex
      neutralParticle->DirectionEnd[0] = annihilationVtx->Direction[0];
      neutralParticle->DirectionEnd[1] = annihilationVtx->Direction[1];
      neutralParticle->DirectionEnd[2] = annihilationVtx->Direction[2];

      // Calculate length
      neutralParticle->Length = sqrt(pow(neutralParticle->PositionEnd[0]-neutralParticle->PositionStart[0],2)+
                                     pow(neutralParticle->PositionEnd[1]-neutralParticle->PositionStart[1],2)+
                                     pow(neutralParticle->PositionEnd[2]-neutralParticle->PositionStart[2],2));

      // Calculate impact parameter
      neutralParticle->ImpactParameter = -999;
      double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
      std::vector<double> parentLineParams;
      pdAnaUtils::ExtrapolateTrack(creationVtx->BeamParticle, parentLineParams, trackFitLength, false);

      if(parentLineParams.size() == 6 && parentLineParams[0] != -999.0){
        TVector3 vertexPos(annihilationVtx->Position[0], annihilationVtx->Position[1], annihilationVtx->Position[2]);
        neutralParticle->ImpactParameter = pdAnaUtils::CalculateImpactParameter(parentLineParams, vertexPos);
      }

      // Calculate neutral particle score
      auto [NPotentialParents, NRecoHitsInVertex] = CalculateNeutralScore(
          neutralParticle, annihilationVtx, creationVtx->BeamParticle, event, particleByUniqueID);

      neutralParticle->AnnihilationVertex->NPotentialParents = NPotentialParents;
      neutralParticle->NRecoHitsInVertex = NRecoHitsInVertex;

      // Fill true equivalent
      neutralParticle->TrueEquivalentNeutralParticle =
          FillTrueEquivalentNeutralParticle(annihilationVtx, creationVtx->BeamParticle);

      // Assign TrueObject
      // Only assign when both true vertex particles have the same true parent (same true ID)
      // AND the parent of the daughters' true parent matches the true object of neutralParticle->Parent
      if(annihilationVtx->TrueEquivalentVertex && annihilationVtx->TrueEquivalentVertex->TrueParticles.size() >= 2){
        AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(annihilationVtx->TrueEquivalentVertex->TrueParticles[0]);
        AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(annihilationVtx->TrueEquivalentVertex->TrueParticles[1]);

        if(trueDaughter1 && trueDaughter2){
          // Both daughters have the same parent (same true ID)
          if(trueDaughter1->ParentID == trueDaughter2->ParentID && trueDaughter1->ParentID > 0){
            // Directly get the true parent particle using the ParentID
            AnaTrueParticlePD* trueParentParticle = pdAnaUtils::GetTrueParticle(&event, trueDaughter1->ParentID);

            if(trueParentParticle){
              // Verify that the parent of the true parent particle matches
              // the true object associated to neutralParticle->Parent
              bool parentMatch = false;

              if(neutralParticle->Parent){
                AnaTrueParticlePD* parentTrueObject = static_cast<AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);

                if(parentTrueObject && trueParentParticle->ParentID > 0){
                  // Get the true grandparent (parent of the daughters' true parent)
                  AnaTrueParticlePD* trueGrandParent = pdAnaUtils::GetTrueParticle(&event, trueParentParticle->ParentID);

                  if(trueGrandParent && trueGrandParent->ID == parentTrueObject->ID){
                    parentMatch = true;
                  }
                }
              }

              if(parentMatch){
                // Assign the true parent to the neutral particle
                neutralParticle->TrueObject = trueParentParticle;
                // Check if the true object has a reconstructed particle associated to it
                if(trueParentParticle->ReconParticles.size() > 0){
                  neutralParticle->RecoParticle = static_cast<AnaParticlePD*>(trueParentParticle->ReconParticles[0]);
                }
              }
              else{
                // Parent match failed - don't assign
                neutralParticle->TrueObject = nullptr;
                neutralParticle->RecoParticle = nullptr;
              }
            }
            else{
              neutralParticle->TrueObject = nullptr;
              neutralParticle->RecoParticle = nullptr;
            }
          }
          else{
            // Daughters have different parents or no parent - don't assign
            neutralParticle->TrueObject = nullptr;
            neutralParticle->RecoParticle = nullptr;
          }
        }
        else{
          neutralParticle->TrueObject = nullptr;
          neutralParticle->RecoParticle = nullptr;
        }
      }
      else{
        neutralParticle->TrueObject = nullptr;
        neutralParticle->RecoParticle = nullptr;
      }

      // Ensure particle momentum
      annihilationVtx->EnsureParticleMomentum();

      // Calculate invariant mass
      Float_t invariantMass = -999;
      if (annihilationVtx->Particles.size() >= 2) {
        AnaParticlePD* particle1 = annihilationVtx->Particles[0];
        AnaParticlePD* particle2 = annihilationVtx->Particles[1];

        if (particle1 && particle2 &&
            particle1->Momentum > 0 && particle2->Momentum > 0 &&
            particle1->Momentum != -999 && particle2->Momentum != -999) {
          const Float_t pionMass = 0.13957;
          invariantMass = anaUtils::ComputeInvariantMass(*particle1, *particle2, pionMass, pionMass);
        }
      }

      neutralParticle->Mass = invariantMass;
      neutralParticle->Momentum = -999;
      neutralParticle->PDG = -999;
      neutralParticle->Lifetime = -999;
      neutralParticle->DecayLength = -999;

      allNeutralParticles.push_back(neutralParticle);
    }
  }

  // Filter out neutral particles whose direction along the beam axis is not positive
  std::vector<AnaNeutralParticlePD*> filteredByDirection;
  for (AnaNeutralParticlePD* neutralParticle : allNeutralParticles) {
    if (!neutralParticle || !neutralParticle->CreationVertex || !neutralParticle->CreationVertex->BeamParticle) {
      delete neutralParticle;
      continue;
    }

    AnaParticlePD* beamParticle = neutralParticle->CreationVertex->BeamParticle;

    // Get beam particle direction (use DirectionEnd for beam)
    TVector3 beamDir(beamParticle->DirectionEnd[0],
                     beamParticle->DirectionEnd[1],
                     beamParticle->DirectionEnd[2]);

    // Get neutral direction
    TVector3 neutralDir(neutralParticle->DirectionStart[0],
                        neutralParticle->DirectionStart[1],
                        neutralParticle->DirectionStart[2]);

    // Normalize directions
    if (beamDir.Mag() > 0 && neutralDir.Mag() > 0) {
      beamDir = beamDir.Unit();
      neutralDir = neutralDir.Unit();

      // Calculate dot product: positive means forward along beam direction
      double alignment = beamDir.Dot(neutralDir);

      // Only keep if alignment is positive (forward along beam)
      if (alignment > 0.0) {
        filteredByDirection.push_back(neutralParticle);
      } else {
        delete neutralParticle;
      }
    } else {
      delete neutralParticle;
    }
  }

  // Normalize scores using alignment with beam particle direction
  NormalizeNeutralParticleScores(filteredByDirection);

  // Filter neutral particles by score: keep only one per annihilation vertex (lowest NeutralScore)
  return FilterNeutralsByScore(filteredByDirection);
}

} // namespace pdNeutralUtils

