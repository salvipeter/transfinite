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
  // TODO
  return Vector3D(0,0,0);
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

Vector3D
RibbonNSided::normal(double s) const {
  double u = inrange(0, s, 1);
  Vector3D n1 = rmf_.eval(0), n2 = rmf_.eval(1);
  return n1 * hermite(0, u) + n2 * hermite(3, u);
}

} // namespace Transfinite
