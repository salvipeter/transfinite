#pragma once

#include <string>

#include "surface-generalized-bezier.hh"
#include "surface-superd.hh"

using namespace Transfinite;

CurveVector readLOP(std::string filename);
TriMesh readOBJ(const std::string &filename);
void writePCP(const std::string &filename, const Point2DVector &uvs, const PointVector &points);
SurfaceGeneralizedBezier loadBezier(const std::string &filename);
void saveBezier(const SurfaceGeneralizedBezier &surf, const std::string &filename);
void writeBezierControlPoints(const SurfaceGeneralizedBezier &surf, const std::string &filename);
std::vector<SurfaceSuperD> loadSuperDModel(const std::string &filename);
