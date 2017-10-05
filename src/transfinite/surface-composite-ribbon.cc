#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-composite-ribbon.hh"
#include "utilities.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceCompositeRibbon::SurfaceCompositeRibbon() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceCompositeRibbon::~SurfaceCompositeRibbon() {
}

Point3D
SurfaceCompositeRibbon::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += compositeRibbon(i, sds[i]) * (blends[i] + blends[prev(i)]);
  return p * 0.5;
}

std::shared_ptr<Ribbon>
SurfaceCompositeRibbon::newRibbon() const {
  return std::make_shared<RibbonType>();
}

Point3D
SurfaceCompositeRibbon::compositeRibbon(size_t i, const Point2D &sd) const {
  double s = sd[0], d = sd[1], s1 = 1.0 - s, d1 = 1.0 - d;
  double Hs = hermite(0, s), Hd = hermite(0, d), Hs1 = 1.0 - Hs;
  return sideInterpolant(i, s, d) * Hd
    + sideInterpolant(prev(i), d1, s) * Hs
    + sideInterpolant(next(i), d, s1) * Hs1
    - cornerCorrection(prev(i), d, s) * Hs * Hd
    - cornerCorrection(i, s1, d) * Hs1 * Hd;
}

} // namespace Transfinite
