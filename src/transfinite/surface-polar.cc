#include "domain-angular.hh"
#include "parameterization-polar.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-polar.hh"
#include "utilities.hh"

namespace Transfinite {

using DomainType = DomainAngular;
using ParamType = ParameterizationPolar;
using RibbonType = RibbonCompatibleWithHandler;

SurfacePolar::SurfacePolar() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfacePolar::~SurfacePolar() {
}

Point3D
SurfacePolar::eval(const Point2D &uv) const {
  Point2DVector pds = param_->mapToRibbons(uv);
  Point3D p(0,0,0);
  double d2sum = 0.0;
  // return polarRibbon(0, pds[0]); // test
  for (size_t i = 0; i < n_; ++i) {
    double d2 = std::pow(pds[i][1], 2);
    p += polarRibbon(i, pds[i]) * d2;
    d2sum += d2;
  }
  return p / d2sum;
}

std::shared_ptr<Ribbon>
SurfacePolar::newRibbon() const {
  return std::make_shared<RibbonType>();
}

Point3D SurfacePolar::polarRibbon(size_t i, const Point2D &pd) const {
  // C0 Bezier displacement version
  double psi = pd[0], d = pd[1];
  d = inrange(0, d, 1);         // avoid -epsilon and 1+epsilon
  size_t i1 = next(i);

  // 1. Linear base interpolant
  // Notations:
  // - uppercase letters are in 3D, lowercase in 2D
  // - 1,2 suffixes indicate the two base sides
  // - P for points, D for (unit) directions, L for lengths
  Point3D P = ribbons_[i]->curve()->eval(1);
  Vector3D D1 = ribbons_[i]->curve()->eval(0) - P;
  Vector3D D2 = ribbons_[i1]->curve()->eval(1) - P;
  double l1 = domain_->edgeLength(i), l2 = domain_->edgeLength(i1);
  double L1 = D1.norm(), L2 = D2.norm();
  D1.normalize(); D2.normalize();
  Vector3D N = D1 ^ D2, D1p = N ^ D1;
  double angle = psi * acos(inrange(-1, D1 * D2, 1));
  Vector3D D = D1 * cos(angle) + D1p * sin(angle);
  Point2D p = domain_->vertices()[i], q = param_->inverse(i, Point2D(psi, 0));
  double l = (q - p).norm();
  double L = l * ((1 - psi) * L1 / l1 + psi * L2 / l2);
  Point3D Q = P + D * L;

  PointVector curve, curve_old = {P, Q};
  curve = curve_old;            // TODO

  return BSCurve(curve).eval(1 - d);
}

// Point3D SurfacePolar::polarRibbon(size_t i, const Point2D &pd) const {
//   double psi = pd[0], d = pd[1];
//   d = inrange(0, d, 1)          // avoid -epsilon and 1+epsilon
//   double d1 = 1.0 - d;
//   return ribbons_[i]->curve()->eval(d) * hermite(0, psi) +
//     ribbons_[i]->crossDerivative(d) * d1 * hermite(1, psi) +
//     ribbons_[next(i)]->crossDerivative(d1) * d1 * hermite(2, psi) +
//     ribbons_[next(i)]->curve()->eval(d1) * hermite(3, psi);
// }

} // namespace Transfinite
