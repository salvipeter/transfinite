#include "domain-angular.hh"
#include "parameterization-parallel.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-side-based.hh"

namespace Transfinite {

using DomainType = DomainAngular;
using ParamType = ParameterizationParallel;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceSideBased::SurfaceSideBased() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceSideBased::~SurfaceSideBased() {
}

Point3D
SurfaceSideBased::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendSideSingular(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += sideInterpolant(i, sds[i][0], sds[i][1]) * blends[i];
  return p;
}

std::shared_ptr<Ribbon>
SurfaceSideBased::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
