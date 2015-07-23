#include "domain-regular.hh"
#include "parameterization-interconnected.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-generalized-coons.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationInterconnected;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceGeneralizedCoons::SurfaceGeneralizedCoons() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceGeneralizedCoons::~SurfaceGeneralizedCoons() {
}

Point3D
SurfaceGeneralizedCoons::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i) {
    double s = sds[i][0], d = sds[i][1], s1 = sds[next(i)][0];
    p += sideInterpolant(i, s, d) * (blends[i] + blends[prev(i)]);
    p -= cornerCorrection(i, 1.0 - s, s1) * blends[i];
  }
  return p;
}

std::shared_ptr<Ribbon>
SurfaceGeneralizedCoons::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
