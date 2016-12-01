#include "domain-regular.hh"
#include "parameterization-bilinear.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-corner-based.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBilinear;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceCornerBased::SurfaceCornerBased() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceCornerBased::~SurfaceCornerBased() {
}

Point3D
SurfaceCornerBased::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += cornerInterpolant(i, sds) * blends[i];
  return p;
}

std::shared_ptr<Ribbon>
SurfaceCornerBased::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
