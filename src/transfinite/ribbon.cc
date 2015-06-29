#include "ribbon.hh"

Ribbon::~Ribbon() {
}

std::shared_ptr<BSCurve>
Ribbon::curve() const {
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
