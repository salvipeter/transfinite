#include "ribbon.hh"

namespace Transfinite {

Ribbon::Ribbon() : multiplier_(1.0), handler_initialized_(false) {
}

Ribbon::~Ribbon() {
}

std::shared_ptr<const Curve>
Ribbon::curve() const {
  return curve_;
}

std::shared_ptr<Curve>
Ribbon::curve() {
  return curve_;
}

void
Ribbon::setCurve(const std::shared_ptr<Curve> &curve) {
  curve_ = curve;
}

void
Ribbon::setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next) {
  prev_ = prev;
  next_ = next;
}

double
Ribbon::multiplier() const {
  return multiplier_;
}

void
Ribbon::setMultiplier(double m) {
  multiplier_ = m;
}

std::optional<Vector3D>
Ribbon::handler() const {
  if (handler_initialized_)
    return handler_;
  return {};
}

void
Ribbon::setHandler(const Vector3D &h) {
  handler_ = h;
  handler_.normalize();
  handler_initialized_ = true;
}

void
Ribbon::overrideNormalFence(const std::shared_ptr<NormalFence> &fence) {
  normal_fence_ = fence;
}

void
Ribbon::reset() {
  multiplier_ = 1.0;
  handler_initialized_ = false;
}

void
Ribbon::update() {
  Vector3D normal;
  VectorVector der;
  rmf_.setCurve(curve_);

  prev_.lock()->curve_->eval(1.0, 1, der);
  normal = der[1];
  curve_->eval(0.0, 1, der);
  normal = normal ^ der[1];
  normal.normalize();
  rmf_.setStart(normal);

  curve_->eval(1.0, 1, der);
  normal = der[1];
  next_.lock()->curve_->eval(0.0, 1, der);
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
  if (normal_fence_)
    return normal_fence_->operator()(s);
  return rmf_.eval(s);
}

} // namespace Transfinite
