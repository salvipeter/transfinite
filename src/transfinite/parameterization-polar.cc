#include "domain.hh"
#include "utilities.hh"
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
  double phi = std::acos(inrange(-1, v1.normalize() * v2.normalize(), 1));
  return Point2D(phi / domain_->angle(i), barycentric(uv)[i]);
}

Point2D
ParameterizationPolar::inverse(size_t i, const Point2D &pd) const {
  if (pd[0] < epsilon)
    return domain_->edgePoint(i, pd[1]);
  else if (pd[0] > 1 - epsilon)
    return domain_->edgePoint(next(i), 1 - pd[1]);

  const size_t precision = 20;  // 2^-20 ~ 1e-6 precision

  Vector2D v1 = domain_->vertices()[prev(i)] - domain_->vertices()[i];
  Vector2D v2 = domain_->vertices()[next(i)] - domain_->vertices()[i];
  v1.normalize(); v2.normalize();
  double phi = pd[0] * domain_->angle(i);
  Vector2D v1p(v1[1], -v1[0]);
  if (v1p * v2 < -epsilon)
    v1p *= -1;
  Point2D p = domain_->vertices()[i], q;
  Vector2D d = v1 * cos(phi) + v1p * sin(phi);

  // Intersect far sides with the sweepline
  for (size_t j = next(next(i)); j != i; j = next(j))
    if (domain_->intersectEdgeWithRay(j, p, d, q))
      break;

  double lo = 0, hi = 1;
  Point2D next;
  for (size_t k = 0; k < precision; ++k) {
    double mid = (lo + hi) / 2;
    next = p * (1 - mid) + q * mid;
    if (mapToRibbon(i, next)[1] > pd[1])
      lo = mid;
    else
      hi = mid;
  }

  return next;
}

} // namespace Transfinite
