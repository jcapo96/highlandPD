#include "neutralKaonTruthTree.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include "BasicUtils.hxx"
#include "TVector3.h"
#include <cmath>

//********************************************************************
void neutralKaonTruthTree::AddNeutralKaonTruthVariables(OutputManager& output, UInt_t nmax){
    AddVarMaxSizeVI(output, k0truepdg, "K0 true PDG", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0truendau, "K0 true number of daughters", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0trueproc, "K0 true process", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0trueendproc, "K0 true end process", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0truestartpos, "K0 true start position", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0trueendpos, "K0 true end position", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0truestartdir, "K0 true start direction", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0trueenddir, "K0 true end direction", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0truestartmom, "K0 true start momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0trueendmom, "K0 true end momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0truelength, "K0 true length", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0truestartenddir, "K0 true start-end direction scalar product", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0truegeneration, "K0 true generation", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0hasrecoobject, "K0 true has reco object", ntruek0, nmax);

    // Vertex reconstruction debugging variables
    AddVarMaxSizeVI(output, k0dau1hasreco, "K0 daughter1 has reco particle", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2hasreco, "K0 daughter2 has reco particle", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parenthasreco, "K0 parent has reco particle", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1validstartpos, "K0 daughter1 reco has valid start position", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1validendpos, "K0 daughter1 reco has valid end position", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2validstartpos, "K0 daughter2 reco has valid start position", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2validendpos, "K0 daughter2 reco has valid end position", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0daucloseenough, "K0 reco daughters are close enough", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0daudistance, "K0 reco daughters distance", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1fitok, "K0 daughter1 track fit succeeded", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2fitok, "K0 daughter2 track fit succeeded", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0vtxpositionfound, "K0 vertex position was found", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0vtxmindistance, "K0 vertex minimum distance between fitted lines", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0parentrecodist, "K0 distance from parent reco end to vertex", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parentwithinradius, "K0 parent within vertex radius", ntruek0, nmax);
}

//********************************************************************
void neutralKaonTruthTree::AddNeutralKaonParentTruthVariables(OutputManager& output, UInt_t nmax){
    AddVarMaxSizeVI(output, k0partruepdg, "K0 parent true PDG", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0partruendau, "K0 parent true number of daughters", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0partrueproc, "K0 parent true process", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0partrueendproc, "K0 parent true end process", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0partruestartpos, "K0 parent true start position", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0partrueendpos, "K0 parent true end position", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0partruestartdir, "K0 parent true start direction", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0partrueenddir, "K0 parent true end direction", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partruestartmom, "K0 parent true start momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partrueendmom, "K0 parent true end momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partruelength, "K0 parent true length", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partruestartenddir, "K0 parent true start-end direction scalar product", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0partruegeneration, "K0 parent true generation", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parhasrecoobject, "K0 parent true has reco object", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parisbeam, "K0 parent true is beam", ntruek0, nmax);
}

//********************************************************************
void neutralKaonTruthTree::AddNeutralKaonDaughter1TruthVariables(OutputManager& output, UInt_t nmax){
    AddVarMaxSizeVI(output, k0dau1truepdg, "K0 daughter1 true PDG", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1truendau, "K0 daughter1 true number of daughters", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1trueproc, "K0 daughter1 true process", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1trueendproc, "K0 daughter1 true end process", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0dau1truestartpos, "K0 daughter1 true start position", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0dau1trueendpos, "K0 daughter1 true end position", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0dau1truestartdir, "K0 daughter1 true start direction", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0dau1trueenddir, "K0 daughter1 true end direction", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1truestartmom, "K0 daughter1 true start momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1trueendmom, "K0 daughter1 true end momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1truelength, "K0 daughter1 true length", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1truestartenddir, "K0 daughter1 true start-end direction scalar product", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1truegeneration, "K0 daughter1 true generation", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1hasrecoobject, "K0 daughter1 true has reco object", ntruek0, nmax);
}

//********************************************************************
void neutralKaonTruthTree::AddNeutralKaonDaughter2TruthVariables(OutputManager& output, UInt_t nmax){
    AddVarMaxSizeVI(output, k0dau2truepdg, "K0 daughter2 true PDG", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2truendau, "K0 daughter2 true number of daughters", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2trueproc, "K0 daughter2 true process", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2trueendproc, "K0 daughter2 true end process", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0dau2truestartpos, "K0 daughter2 true start position", ntruek0, nmax);
    AddVarMaxSize4MF(output, k0dau2trueendpos, "K0 daughter2 true end position", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0dau2truestartdir, "K0 daughter2 true start direction", ntruek0, nmax);
    AddVarMaxSize3MF(output, k0dau2trueenddir, "K0 daughter2 true end direction", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2truestartmom, "K0 daughter2 true start momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2trueendmom, "K0 daughter2 true end momentum", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2truelength, "K0 daughter2 true length", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2truestartenddir, "K0 daughter2 true start-end direction scalar product", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2truegeneration, "K0 daughter2 true generation", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2hasrecoobject, "K0 daughter2 true has reco object", ntruek0, nmax);
}

//********************************************************************
void neutralKaonTruthTree::FillNeutralKaonTruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject){
    output.FillVectorVar(k0truepdg, part.PDG);
    output.FillVectorVar(k0truendau, static_cast<Int_t>(part.Daughters.size()));
    output.FillVectorVar(k0trueproc, part.ProcessStart);
    output.FillVectorVar(k0trueendproc, part.ProcessEnd);
    output.FillMatrixVarFromArray(k0truestartpos, part.Position, 4);
    output.FillMatrixVarFromArray(k0trueendpos, part.PositionEnd, 4);
    output.FillMatrixVarFromArray(k0truestartdir, part.Direction, 3);
    output.FillMatrixVarFromArray(k0trueenddir, part.DirectionEnd, 3);
    output.FillVectorVar(k0truestartmom, part.Momentum);
    output.FillVectorVar(k0trueendmom, part.MomentumEnd);

    //check if the reco object exists
    output.FillVectorVar(k0hasrecoobject, hasRecoObject ? 1 : 0);

    // Calculate length from positions
    Float_t length = sqrt(pow(part.PositionEnd[0] - part.Position[0], 2) +
                         pow(part.PositionEnd[1] - part.Position[1], 2) +
                         pow(part.PositionEnd[2] - part.Position[2], 2));
    output.FillVectorVar(k0truelength, length);

    // Calculate dot product of start and end directions
    Float_t dotProduct = part.Direction[0] * part.DirectionEnd[0] +
                        part.Direction[1] * part.DirectionEnd[1] +
                        part.Direction[2] * part.DirectionEnd[2];
    output.FillVectorVar(k0truestartenddir, dotProduct);
    output.FillVectorVar(k0truegeneration, part.Generation);
}

//********************************************************************
void neutralKaonTruthTree::FillVertexReconstructionDebugVariables(OutputManager& output, const AnaTrueParticlePD& part,
                                                                  AnaParticlePD* daughter1Reco, AnaParticlePD* daughter2Reco,
                                                                  AnaParticlePD* parentReco,
                                                                  double maxDaughterDistance, double trackFitLength, double vertexRadius){

    // Initialize all debug variables to 0 or -999
    output.FillVectorVar(k0dau1hasreco, daughter1Reco != nullptr ? 1 : 0);
    output.FillVectorVar(k0dau2hasreco, daughter2Reco != nullptr ? 1 : 0);
    output.FillVectorVar(k0parenthasreco, parentReco != nullptr ? 1 : 0);

    if(!daughter1Reco || !daughter2Reco){
        // If either daughter doesn't have a reco particle, fill remaining variables with defaults
        output.FillVectorVar(k0dau1validstartpos, 0);
        output.FillVectorVar(k0dau1validendpos, 0);
        output.FillVectorVar(k0dau2validstartpos, 0);
        output.FillVectorVar(k0dau2validendpos, 0);
        output.FillVectorVar(k0daucloseenough, 0);
        output.FillVectorVar(k0daudistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0dau1fitok, 0);
        output.FillVectorVar(k0dau2fitok, 0);
        output.FillVectorVar(k0vtxpositionfound, 0);
        output.FillVectorVar(k0vtxmindistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentrecodist, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentwithinradius, 0);
        return;
    }

    // Check if daughter1 has valid start position
    bool dau1ValidStart = (daughter1Reco->PositionStart[0] > -900 &&
                          daughter1Reco->PositionStart[1] > -900 &&
                          daughter1Reco->PositionStart[2] > -900);
    output.FillVectorVar(k0dau1validstartpos, dau1ValidStart ? 1 : 0);

    // Check if daughter1 has valid end position
    bool dau1ValidEnd = (daughter1Reco->PositionEnd[0] > -900 &&
                        daughter1Reco->PositionEnd[1] > -900 &&
                        daughter1Reco->PositionEnd[2] > -900);
    output.FillVectorVar(k0dau1validendpos, dau1ValidEnd ? 1 : 0);

    // Check if daughter2 has valid start position
    bool dau2ValidStart = (daughter2Reco->PositionStart[0] > -900 &&
                          daughter2Reco->PositionStart[1] > -900 &&
                          daughter2Reco->PositionStart[2] > -900);
    output.FillVectorVar(k0dau2validstartpos, dau2ValidStart ? 1 : 0);

    // Check if daughter2 has valid end position
    bool dau2ValidEnd = (daughter2Reco->PositionEnd[0] > -900 &&
                        daughter2Reco->PositionEnd[1] > -900 &&
                        daughter2Reco->PositionEnd[2] > -900);
    output.FillVectorVar(k0dau2validendpos, dau2ValidEnd ? 1 : 0);

    if(!dau1ValidStart || !dau1ValidEnd || !dau2ValidStart || !dau2ValidEnd){
        output.FillVectorVar(k0daucloseenough, 0);
        output.FillVectorVar(k0daudistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0dau1fitok, 0);
        output.FillVectorVar(k0dau2fitok, 0);
        output.FillVectorVar(k0vtxpositionfound, 0);
        output.FillVectorVar(k0vtxmindistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentrecodist, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentwithinradius, 0);
        return;
    }

    // Check if daughters are close enough
    TVector3 pos1 = pdAnaUtils::DefinePosition(daughter1Reco);
    TVector3 pos2 = pdAnaUtils::DefinePosition(daughter2Reco);

    Float_t pos1_array[3] = {static_cast<Float_t>(pos1.X()), static_cast<Float_t>(pos1.Y()), static_cast<Float_t>(pos1.Z())};
    Float_t pos2_array[3] = {static_cast<Float_t>(pos2.X()), static_cast<Float_t>(pos2.Y()), static_cast<Float_t>(pos2.Z())};
    Float_t distance = static_cast<Float_t>(sqrt(anaUtils::GetSeparationSquared(pos1_array, pos2_array)));

    output.FillVectorVar(k0daudistance, static_cast<Float_t>(distance));
    bool closeEnough = (distance <= static_cast<Float_t>(maxDaughterDistance));
    output.FillVectorVar(k0daucloseenough, closeEnough ? 1 : 0);

    if(!closeEnough){
        output.FillVectorVar(k0dau1fitok, 0);
        output.FillVectorVar(k0dau2fitok, 0);
        output.FillVectorVar(k0vtxpositionfound, 0);
        output.FillVectorVar(k0vtxmindistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentrecodist, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentwithinradius, 0);
        return;
    }

    // Try to fit tracks
    std::vector<double> line1Params, line2Params;
    pdAnaUtils::ExtrapolateTrack(daughter1Reco, line1Params, trackFitLength, true);
    pdAnaUtils::ExtrapolateTrack(daughter2Reco, line2Params, trackFitLength, true);

    bool dau1FitOk = (line1Params[0] != -999.0);
    bool dau2FitOk = (line2Params[0] != -999.0);
    output.FillVectorVar(k0dau1fitok, dau1FitOk ? 1 : 0);
    output.FillVectorVar(k0dau2fitok, dau2FitOk ? 1 : 0);

    if(!dau1FitOk || !dau2FitOk){
        output.FillVectorVar(k0vtxpositionfound, 0);
        output.FillVectorVar(k0vtxmindistance, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentrecodist, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentwithinradius, 0);
        return;
    }

    // Find closest points between fitted lines
    TVector3 closestPoint1, closestPoint2;
    double minDistance = pdAnaUtils::FindClosestPointsBetweenLines(line1Params, line2Params,
                                                                   closestPoint1, closestPoint2);
    output.FillVectorVar(k0vtxmindistance, static_cast<Float_t>(minDistance));

    // Calculate vertex position
    TVector3 vertexPosition = 0.5 * (closestPoint1 + closestPoint2);
    bool vtxFound = (vertexPosition.X() > -900 && vertexPosition.Y() > -900 && vertexPosition.Z() > -900);
    output.FillVectorVar(k0vtxpositionfound, vtxFound ? 1 : 0);

    // Check distance from parent end position to vertex
    if(parentReco && vtxFound){
        TVector3 parentEnd(parentReco->PositionEnd[0], parentReco->PositionEnd[1], parentReco->PositionEnd[2]);
        Float_t parentDist = static_cast<Float_t>((vertexPosition - parentEnd).Mag());
        output.FillVectorVar(k0parentrecodist, parentDist);

        bool withinRadius = (parentDist <= static_cast<Float_t>(vertexRadius));
        output.FillVectorVar(k0parentwithinradius, withinRadius ? 1 : 0);
    }
    else{
        std::cout << "DEBUG: No parent or vertex found" << std::endl;
        std::cout << "DEBUG: Parent ID: " << part.ID << std::endl;
        std::cout << "DEBUG: Vertex position: " << vertexPosition.X() << ", " << vertexPosition.Y() << ", " << vertexPosition.Z() << std::endl;
        output.FillVectorVar(k0parentrecodist, static_cast<Float_t>(-999.0));
        output.FillVectorVar(k0parentwithinradius, 0);
    }
}

//********************************************************************
void neutralKaonTruthTree::FillNeutralKaonParentTruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject){
    output.FillVectorVar(k0partruepdg, part.PDG);
    output.FillVectorVar(k0partruendau, static_cast<Int_t>(part.Daughters.size()));
    output.FillVectorVar(k0partrueproc, part.ProcessStart);
    output.FillVectorVar(k0partrueendproc, part.ProcessEnd);
    output.FillMatrixVarFromArray(k0partruestartpos, part.Position, 4);
    output.FillMatrixVarFromArray(k0partrueendpos, part.PositionEnd, 4);
    output.FillMatrixVarFromArray(k0partruestartdir, part.Direction, 3);
    output.FillMatrixVarFromArray(k0partrueenddir, part.DirectionEnd, 3);
    output.FillVectorVar(k0partruestartmom, part.Momentum);
    output.FillVectorVar(k0partrueendmom, part.MomentumEnd);
    // output.FillVectorVar(k0parisbeam, part.IsBeamPart ? 1 : 0);

    // Calculate length from positions
    Float_t length = sqrt(pow(part.PositionEnd[0] - part.Position[0], 2) +
                         pow(part.PositionEnd[1] - part.Position[1], 2) +
                         pow(part.PositionEnd[2] - part.Position[2], 2));
    output.FillVectorVar(k0partruelength, length);

    // Calculate dot product of start and end directions
    Float_t dotProduct = part.Direction[0] * part.DirectionEnd[0] +
                        part.Direction[1] * part.DirectionEnd[1] +
                        part.Direction[2] * part.DirectionEnd[2];
    output.FillVectorVar(k0partruestartenddir, dotProduct);
    output.FillVectorVar(k0partruegeneration, part.Generation);

    // Check if the reco object exists
    output.FillVectorVar(k0parhasrecoobject, hasRecoObject ? 1 : 0);
}

//********************************************************************
void neutralKaonTruthTree::FillNeutralKaonDaughter1TruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject){
    output.FillVectorVar(k0dau1hasrecoobject, hasRecoObject ? 1 : 0);
    output.FillVectorVar(k0dau1truepdg, part.PDG);
    output.FillVectorVar(k0dau1truendau, static_cast<Int_t>(part.Daughters.size()));
    output.FillVectorVar(k0dau1trueproc, part.ProcessStart);
    output.FillVectorVar(k0dau1trueendproc, part.ProcessEnd);
    output.FillMatrixVarFromArray(k0dau1truestartpos, part.Position, 4);
    output.FillMatrixVarFromArray(k0dau1trueendpos, part.PositionEnd, 4);
    output.FillMatrixVarFromArray(k0dau1truestartdir, part.Direction, 3);
    output.FillMatrixVarFromArray(k0dau1trueenddir, part.DirectionEnd, 3);
    output.FillVectorVar(k0dau1truestartmom, part.Momentum);
    output.FillVectorVar(k0dau1trueendmom, part.MomentumEnd);

    // Calculate length from positions
    Float_t length = sqrt(pow(part.PositionEnd[0] - part.Position[0], 2) +
                         pow(part.PositionEnd[1] - part.Position[1], 2) +
                         pow(part.PositionEnd[2] - part.Position[2], 2));
    output.FillVectorVar(k0dau1truelength, length);

    // Calculate dot product of start and end directions
    Float_t dotProduct = part.Direction[0] * part.DirectionEnd[0] +
                        part.Direction[1] * part.DirectionEnd[1] +
                        part.Direction[2] * part.DirectionEnd[2];
    output.FillVectorVar(k0dau1truestartenddir, dotProduct);
    output.FillVectorVar(k0dau1truegeneration, part.Generation);
}

//********************************************************************
void neutralKaonTruthTree::FillNeutralKaonDaughter2TruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject){
    output.FillVectorVar(k0dau2hasrecoobject, hasRecoObject ? 1 : 0);
    output.FillVectorVar(k0dau2truepdg, part.PDG);
    output.FillVectorVar(k0dau2truendau, static_cast<Int_t>(part.Daughters.size()));
    output.FillVectorVar(k0dau2trueproc, part.ProcessStart);
    output.FillVectorVar(k0dau2trueendproc, part.ProcessEnd);
    output.FillMatrixVarFromArray(k0dau2truestartpos, part.Position, 4);
    output.FillMatrixVarFromArray(k0dau2trueendpos, part.PositionEnd, 4);
    output.FillMatrixVarFromArray(k0dau2truestartdir, part.Direction, 3);
    output.FillMatrixVarFromArray(k0dau2trueenddir, part.DirectionEnd, 3);
    output.FillVectorVar(k0dau2truestartmom, part.Momentum);
    output.FillVectorVar(k0dau2trueendmom, part.MomentumEnd);

    // Calculate length from positions
    Float_t length = sqrt(pow(part.PositionEnd[0] - part.Position[0], 2) +
                         pow(part.PositionEnd[1] - part.Position[1], 2) +
                         pow(part.PositionEnd[2] - part.Position[2], 2));
    output.FillVectorVar(k0dau2truelength, length);

    // Calculate dot product of start and end directions
    Float_t dotProduct = part.Direction[0] * part.DirectionEnd[0] +
                        part.Direction[1] * part.DirectionEnd[1] +
                        part.Direction[2] * part.DirectionEnd[2];
    output.FillVectorVar(k0dau2truestartenddir, dotProduct);
    output.FillVectorVar(k0dau2truegeneration, part.Generation);
}