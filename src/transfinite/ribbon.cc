#include "ribbon.hh"

namespace Transfinite {

Ribbon::~Ribbon() {
}

std::shared_ptr<const BSCurve>
Ribbon::curve() const {
  return curve_;
}

std::shared_ptr<BSCurve>
Ribbon::curve() {
  return curve_;
}

void
Ribbon::setCurve(const std::shared_ptr<BSCurve> &curve) {
  curve_ = curve;
}

void
Ribbon::setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next) {
  prev_ = prev;
  next_ = next;
}

void
Ribbon::update() {
  Vector3D normal;
  VectorVector der;
  rmf_.setCurve(curve_);

  prev_->curve_->eval(1.0, 1, der);
  normal = der[1];
  curve_->eval(0.0, 1, der);
  normal = normal ^ der[1];
  normal.normalize();
  rmf_.setStart(normal);

  curve_->eval(1.0, 1, der);
  normal = der[1];
  next_->curve_->eval(0.0, 1, der);
  normal = normal ^ der[1];
  normal.normalize();
  rmf_.setEnd(normal);

  rmf_.update();
}

Point3D
Ribbon::eval(const Point2D &sd) const {
  return curve_->eval(sd[0]) + crossDerivative(sd[0]) * sd[1];
}

Vector3D
Ribbon::normal(double s) const {
  return rmf_.eval(s);
}

} // namespace Transfinite
