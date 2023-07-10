#include <iostream>
#include "ribbon-coons.hh"
#include "utilities.hh"

namespace Transfinite {

RibbonCoons::~RibbonCoons() {
}

void
RibbonCoons::update() {
  left_ = prev_.lock()->curve();
  right_ = next_.lock()->curve();
  auto left2 = dynamic_cast<RibbonCoons *>(prev_.lock().get())->prev_.lock()->curve();
  auto right2 = dynamic_cast<RibbonCoons *>(next_.lock().get())->next_.lock()->curve();
  auto p1 = left_->eval(0.0);
  auto q1 = right_->eval(1.0);
  if (right2 == left_)          // 3-sided
    top_ = std::make_shared<BSplineCurve>(BSCurve(PointVector{p1})); // 0-degree Bezier curve
  else {
    VectorVector der;
    left2->eval(1.0, 1, der);
    auto p2 = p1 - der[1] / 3.0;
    right2->eval(0.0, 1, der);
    auto q2 = q1 + der[1] / 3.0;
    top_ = std::make_shared<BSplineCurve>(BSCurve(PointVector{q1, q2, p2, p1}));
  }
  bl_ = curve_->eval(0.0);
  br_ = curve_->eval(1.0);
  tl_ = left_->eval(0.0);
  tr_ = right_->eval(1.0);
  Ribbon::update();
}

Vector3D
RibbonCoons::crossDerivative(double) const {
  // TODO
  return {0, 0, 0};
}

Point3D
RibbonCoons::eval(const Point2D &sd) const {
  auto s = inrange(0, sd[0], 1), d = inrange(0, sd[1], 1);
  auto s1 = inrange(0, 1 - s, 1), d1 = inrange(0, 1 - d, 1);
  auto p1 = curve_->eval(s) * d1 + top_->eval(s1) * d;
  auto p2 = left_->eval(d1) * s1 + right_->eval(d) * s;
  auto p12 = (bl_ * s1 + br_ * s) * d1 + (tl_ * s1 + tr_ * s) * d;
  return p1 + p2 - p12;
}

} // namespace Transfinite
