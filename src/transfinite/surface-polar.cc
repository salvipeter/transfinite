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
  double psi = pd[0], d = pd[1];
  d = std::min(std::max(d, 0.0), 1.0); // avoid -epsilon and 1+epsilon
  double d1 = 1.0 - d;
  return ribbons_[i]->curve()->eval(d) * hermite(0, psi) +
    ribbons_[i]->crossDerivative(d) * d1 * hermite(1, psi) +
    ribbons_[next(i)]->crossDerivative(d1) * d1 * hermite(2, psi) +
    ribbons_[next(i)]->curve()->eval(d1) * hermite(3, psi);
}

} // namespace Transfinite
