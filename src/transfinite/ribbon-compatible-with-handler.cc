#include "ribbon-compatible-with-handler.hh"

namespace Transfinite {

RibbonCompatibleWithHandler::~RibbonCompatibleWithHandler() {
}

void
RibbonCompatibleWithHandler::update() {
  RibbonCompatible::update();
  Vector3D n = normal(0.5);
  if (!handler_initialized_) {
    handler_ = prev_tangent_ / prev_tangent_.norm() + next_tangent_ / next_tangent_.norm();
    handler_ = handler_ - n * (handler_ * n);
    handler_.normalize();
  }
  central_ = handler_ * (prev_tangent_.norm() + next_tangent_.norm()) / 2.0 * multiplier_;
}

Vector3D
RibbonCompatibleWithHandler::crossDerivative(double s) const {
  Vector3D n = normal(s);
  Vector3D pt = prev_tangent_ - n * (prev_tangent_ * n);
  Vector3D ch = central_ - n * (central_ * n);
  Vector3D nt = next_tangent_ - n * (next_tangent_ * n);
  return pt * 2.0 * (s - 1.0) * (s - 0.5)
    + ch * -4.0 * s * (s - 1.0)
    + nt * 2.0 * s * (s - 0.5);
}

} // namespace Transfinite
