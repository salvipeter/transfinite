#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-dummy.hh"
#include "surface-generalized-bezier.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonDummy;

SurfaceGeneralizedBezier::SurfaceGeneralizedBezier() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceGeneralizedBezier::~SurfaceGeneralizedBezier() {
}

Point3D
SurfaceGeneralizedBezier::eval(const Point2D &uv) const {
  // TODO
  return Point3D(0,0,0);
}

std::shared_ptr<Ribbon>
SurfaceGeneralizedBezier::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
