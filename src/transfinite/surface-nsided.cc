#include "domain-angular.hh"
#include "parameterization-perp-polar.hh"
#include "ribbon-nsided.hh"
#include "surface-nsided.hh"

namespace Transfinite {

using DomainType = DomainAngular;
using ParamType = ParameterizationPerpPolar;
using RibbonType = RibbonNSided;

SurfaceNSided::SurfaceNSided() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceNSided::~SurfaceNSided() {
}

Point3D
SurfaceNSided::eval(const Point2D &uv) const {
  // Test only one ribbon for now
  size_t ribbon_index = 3;
  Point2D sd = param_->mapToRibbon(ribbon_index, uv);
  return ribbons_[ribbon_index]->eval(sd);
}

std::shared_ptr<Ribbon>
SurfaceNSided::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
