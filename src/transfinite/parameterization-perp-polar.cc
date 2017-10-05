#include <cmath>

#include "domain.hh"
#include "parameterization-perp-polar.hh"
#include "utilities.hh"

namespace Transfinite {

ParameterizationPerpPolar::~ParameterizationPerpPolar() {
}

Point2D
ParameterizationPerpPolar::mapToRibbon(size_t i, const Point2D &uv) const {
  Point2DVector const &vertices = domain_->vertices();

  const Point2D &base = vertices[prev(i)];

  Point2D p = domain_->toLocal(i, uv - base);
  if (p[0] >= 0 && p[0] <= 1)
    return p;

  Point2D origin = base;
  Vector2D vbase = domain_->fromLocal(i, Vector2D(0, 1));
  if (p[0] > 1) {
    origin = vertices[i];
    vbase *= -1;
  }

  Vector2D v = uv - origin;
  double d = v.norm() / domain_->edgeLength(i);
  double phi = std::acos(inrange(-1, vbase.normalize() * v.normalize(), 1));

  if (p[0] < 0)
    return Point2D(-phi, d);
  return Point2D(1+phi, d);
}

} // namespace Transfinite
