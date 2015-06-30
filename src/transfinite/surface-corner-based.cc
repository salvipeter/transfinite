#include "domain-regular.hh"
#include "parameterization-bilinear.hh"
#include "ribbon-compatible.hh"
#include "surface-corner-based.hh"

using DomainType = DomainRegular;
using ParamType = ParameterizationBilinear;
using RibbonType = RibbonCompatible;

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

Point3D
SurfaceCornerBased::cornerInterpolant(size_t i, const Point2DVector &sds) const {
  double si = sds[i][0], si1 = sds[next(i)][0];
  return sideInterpolant(i, si, si1) + sideInterpolant(next(i), si1, 1.0 - si)
    - cornerCorrection(i, 1.0 - si, si1);
}
