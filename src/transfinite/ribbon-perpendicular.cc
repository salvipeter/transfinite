#include "ribbon-perpendicular.hh"

namespace Transfinite {

RibbonPerpendicular::~RibbonPerpendicular() {
}

void
RibbonPerpendicular::update() {
  Ribbon::update();

  VectorVector der;
  curve_->eval(0.0, 1, der);
  auto t0 = der[1].normalize();
  curve_->eval(1.0, 1, der);
  auto t1 = der[1].normalize();
  prev_.lock()->curve()->eval(1.0, 1, der);
  auto pt = -der[1];
  next_.lock()->curve()->eval(0.0, 1, der);
  auto nt = der[1];
  auto n0 = normal(0.0);    auto n1 = normal(1.0);
  auto b0 = n0 ^ t0;        auto b1 = n1 ^ t1;
  prev_norm_ = pt.norm();   next_norm_ = nt.norm();
  pt /= prev_norm_;         nt /= next_norm_;
  prev_alpha_ = pt * b0;    prev_beta_ = pt * t0;
  next_alpha_ = nt * b1;    next_beta_ = nt * t1;
}

Vector3D
RibbonPerpendicular::crossDerivative(double s) const {
  VectorVector der;
  curve_->eval(s, 1, der);
  auto t = der[1].normalize();
  auto n = normal(s);
  auto size = prev_norm_ * (1.0 - s) + next_norm_ * s;
  auto alpha = prev_alpha_ * 2.0 * (s - 1.0) * (s - 0.5) + next_alpha_ * 2.0 * s * (s - 0.5)
    - 4.0 * s * (s - 1.0);
  auto beta = prev_beta_ * 2.0 * (s - 1.0) * (s - 0.5) + next_beta_ * 2.0 * s * (s - 0.5);
    
  return ((n ^ t) * alpha + t * beta) * size;
}

} // namespace Transfinite
