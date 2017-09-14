#include "domain.hh"
#include "parameterization-polar.hh"

namespace Transfinite {

ParameterizationPolar::ParameterizationPolar() :
  ParameterizationBarycentric(BarycentricType::MEAN_VALUE) {
}

ParameterizationPolar::~ParameterizationPolar() {
}

Point2D
ParameterizationPolar::mapToRibbon(size_t i, const Point2D &uv) const {
  Vector2D v1 = domain_->vertices()[prev(i)] - domain_->vertices()[i];
  Vector2D v2 = uv - domain_->vertices()[i];
  if (v2.norm() < epsilon)
    return Point2D(0.0, 1.0);
  double phi = std::acos(std::min(std::max(v1.normalize() * v2.normalize(), -1.0), 1.0));
  return Point2D(phi / domain_->angle(i), barycentric(uv)[i]);
}

} // namespace Transfinite
