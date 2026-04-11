#include "pdCreationUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include <algorithm>
#include <cmath>

namespace pdCreationUtils {

namespace {

bool HasValid3DPoint(const TVector3& point) {
	return std::isfinite(point.X()) && std::isfinite(point.Y()) && std::isfinite(point.Z()) &&
				 point.X() > -900.0 && point.Y() > -900.0 && point.Z() > -900.0;
}

bool HasValidParticlePoint(const Float_t point[4]) {
	if (!point) return false;
	return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]) &&
				 point[0] > -900.0 && point[1] > -900.0 && point[2] > -900.0;
}

} // namespace

TVector3 ProjectPointOntoLine(const TVector3& point,
															const TVector3& linePoint,
															const TVector3& lineDirection) {
	if (!HasValid3DPoint(point) || !HasValid3DPoint(linePoint)) {
		return TVector3(-999.0, -999.0, -999.0);
	}

	TVector3 dir = lineDirection;
	if (dir.Mag2() <= 1e-10 || !std::isfinite(dir.X()) || !std::isfinite(dir.Y()) || !std::isfinite(dir.Z())) {
		return TVector3(-999.0, -999.0, -999.0);
	}

	dir = dir.Unit();
	const TVector3 delta = point - linePoint;
	const double t = delta.Dot(dir);
	return linePoint + t * dir;
}

bool ProjectBeamTailOntoStartDirection(const AnaParticlePD* beamParticle,
									   double fitDistanceFromEndCm,
																			 TVector3& correctedEndPosition,
																			 std::vector<TVector3>* rawTailPoints,
																			 std::vector<TVector3>* projectedTailPoints) {
	correctedEndPosition.SetXYZ(-999.0, -999.0, -999.0);
	if (rawTailPoints) rawTailPoints->clear();
	if (projectedTailPoints) projectedTailPoints->clear();
	if (!beamParticle) return false;

	double fitLengthCm = 25.0;
	double fitDistanceFromEndForFitCm = 10.0;
	if (ND::params().HasParameter("neutralKaonAnalysis.TrackFitCreationLength")) {
		fitLengthCm = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitCreationLength");
	}
	if (ND::params().HasParameter("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd")) {
		fitDistanceFromEndForFitCm = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd");
	}

	std::vector<double> fitParams;
	pdAnaUtils::ExtrapolateTrack(const_cast<AnaParticlePD*>(beamParticle),
															 fitParams,
															 fitLengthCm,
															 false,
															   fitDistanceFromEndForFitCm);

	TVector3 linePoint(-999.0, -999.0, -999.0);
	TVector3 lineDirection(-999.0, -999.0, -999.0);
	const bool hasValidFit =
			(fitParams.size() >= 6 && std::isfinite(fitParams[0]) && std::isfinite(fitParams[1]) &&
			 std::isfinite(fitParams[2]) && std::isfinite(fitParams[3]) && std::isfinite(fitParams[4]) &&
			 std::isfinite(fitParams[5]) && fitParams[0] > -900.0 && fitParams[1] > -900.0 && fitParams[2] > -900.0);

	if (hasValidFit) {
		linePoint.SetXYZ(fitParams[0], fitParams[1], fitParams[2]);
		lineDirection.SetXYZ(fitParams[3], fitParams[4], fitParams[5]);
	} else {
		linePoint.SetXYZ(beamParticle->PositionEnd[0], beamParticle->PositionEnd[1], beamParticle->PositionEnd[2]);
		lineDirection.SetXYZ(-beamParticle->DirectionEnd[0], -beamParticle->DirectionEnd[1], -beamParticle->DirectionEnd[2]);
	}

	if (!HasValid3DPoint(linePoint) || lineDirection.Mag2() <= 1e-10) {
		return false;
	}
	lineDirection = lineDirection.Unit();

	const double effectiveTailLength = std::max(0.0, fitDistanceFromEndCm);

	std::vector<TVector3> tailRaw;
	tailRaw.reserve(beamParticle->TrjPoints.size());

	TVector3 prevPoint;
	bool hasPrevPoint = false;
	double cumulativeLength = 0.0;

	for (int i = static_cast<int>(beamParticle->TrjPoints.size()) - 1; i >= 0; --i) {
		const TVector3& pos = beamParticle->TrjPoints[i].Position;
		if (!HasValid3DPoint(pos)) continue;

		if (hasPrevPoint) {
			cumulativeLength += (prevPoint - pos).Mag();
		}

		if (tailRaw.empty() || cumulativeLength <= effectiveTailLength + 1e-6) {
			tailRaw.push_back(pos);
			prevPoint = pos;
			hasPrevPoint = true;
			continue;
		}

		break;
	}

	if (tailRaw.empty() && HasValidParticlePoint(beamParticle->PositionEnd)) {
		tailRaw.emplace_back(beamParticle->PositionEnd[0], beamParticle->PositionEnd[1], beamParticle->PositionEnd[2]);
	}

	if (tailRaw.empty()) return false;

	std::reverse(tailRaw.begin(), tailRaw.end());
	std::vector<TVector3> tailProjected;
	tailProjected.reserve(tailRaw.size());
	for (const TVector3& rawPoint : tailRaw) {
		const TVector3 projected = ProjectPointOntoLine(rawPoint, linePoint, lineDirection);
		if (!HasValid3DPoint(projected)) continue;
		tailProjected.push_back(projected);
	}

	if (tailProjected.empty()) return false;

	correctedEndPosition = tailProjected.back();

	if (rawTailPoints) *rawTailPoints = tailRaw;
	if (projectedTailPoints) *projectedTailPoints = tailProjected;

	return true;
}
}
