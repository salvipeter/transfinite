#pragma once

#include "surface-generalized-bezier.hh"

using namespace Transfinite;

SurfaceGeneralizedBezier elevateDegree(const SurfaceGeneralizedBezier &surf);
Point2DVector parameterizePoints(const Surface &surf, const PointVector &points);
SurfaceGeneralizedBezier fitWithOriginal(const SurfaceGeneralizedBezier &original,
                                         const PointVector &points,
                                         const Point2DVector &params,
                                         double smoothing = 0, size_t fixed_rows = 2);
