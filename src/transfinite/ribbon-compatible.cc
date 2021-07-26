#include "ribbon-compatible.hh"

namespace Transfinite {

RibbonCompatible::~RibbonCompatible() {
}

void
RibbonCompatible::update() {
  VectorVector der;
  prev_.lock()->curve()->eval(1.0, 1, der);
  prev_tangent_ = -der[1];
  next_.lock()->curve()->eval(0.0, 1, der);
  next_tangent_ = der[1];

  Ribbon::update();
}

Vector3D
RibbonCompatible::crossDerivative(double s) const {
  Vector3D n = normal(s);
  Vector3D pt = prev_tangent_ - n * (prev_tangent_ * n);
  Vector3D nt = next_tangent_ - n * (next_tangent_ * n);
  return pt * (1.0 - s) + nt * s;
}

} // namespace Transfinite
