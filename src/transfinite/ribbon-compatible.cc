#include "ribbon-compatible.hh"

RibbonCompatible::~RibbonCompatible() {
}

void
RibbonCompatible::update() {
  VectorVector der;
  prev_->curve()->eval(1.0, 1, der);
  prev_tangent_ = -der[1];
  next_->curve()->eval(0.0, 1, der);
  next_tangent_ = der[1];

  Ribbon::update();
}

Point3D
RibbonCompatible::eval(const Point2D &sd) const {
  Vector3D n = rmf_.eval(sd[0]);
  Vector3D pt = prev_tangent_ - n * (prev_tangent_ * n);
  Vector3D nt = next_tangent_ - n * (next_tangent_ * n);
  return curve_->eval(sd[0]) + (pt * (1.0 - sd[0]) + nt * sd[0]) * sd[1];
}
