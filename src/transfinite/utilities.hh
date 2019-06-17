#pragma once

#include "geometry.hh"

namespace Transfinite {

using namespace Geometry;

Point3D affineCombine(const Point3D &p, double x, const Point3D &q);

double inrange(double min, double x, double max);

double hermite(int i, double t);

void bernstein(size_t n, double u, DoubleVector &coeff);
double bernstein(size_t i, size_t n, double u);
void bezierElevate(PointVector &cpts);

} // namespace Transfinite
