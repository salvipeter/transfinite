#pragma once

#include "geometry.hh"

namespace Transfinite {

using namespace Geometry;

double inrange(double min, double x, double max);

double blendHermite(double x);
double hermite(int i, double t);

void bernstein(size_t n, double u, DoubleVector &coeff);
double bernstein(size_t i, size_t n, double u);
void bezierElevate(PointVector &cpts);

} // namespace Transfinite
