#include <cmath>

#include "ribbon-nsided.hh"
#include "utilities.hh"

namespace Transfinite {

RibbonNSided::~RibbonNSided() {
}

void
RibbonNSided::update() {
  base_length_ = curve_->arcLength(0, 1);

  Ribbon::update();
}

Vector3D
RibbonNSided::crossDerivative(double s) const {
  Vector3D n = normal(s);
  double u = inrange(0, s, 1);
  VectorVector der;
  curve_->eval(u, 1, der);
  Vector3D &d = der[1].normalize();
  if (s == u)
    return n ^ d;

  // Polar lines
  double phi = s < 0 ? -s : s - 1;
  if (s < 0)
    d *= -1;
  Vector3D b = d ^ n;
  return b * std::cos(phi) + d * std::sin(phi);
}

Point3D
RibbonNSided::eval(const Point2D &sd) const {
  // s ranges from -pi to 1+pi, where values <0 and >1 are meant as angles
  // d is the ratio of the sweep in the domain to the length of the base side
  // d is negative for points "behind" the ribbon
  double u = inrange(0, sd[0], 1);
  Point3D p = curve_->eval(u);
  double length = sd[1] * base_length_;
  return p + crossDerivative(sd[0]) * length;
}

} // namespace Transfinite
