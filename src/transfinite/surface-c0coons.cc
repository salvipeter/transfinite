#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-coons.hh"
#include "surface-c0coons.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonCoons;

SurfaceC0Coons::SurfaceC0Coons() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceC0Coons::~SurfaceC0Coons() {
}

Point3D
SurfaceC0Coons::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += ribbons_[i]->eval(sds[i]) * (1 - sds[i][1]) / 2;
  return p;
}

std::shared_ptr<Ribbon>
SurfaceC0Coons::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
