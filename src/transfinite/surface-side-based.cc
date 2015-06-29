#include "domain-regular.hh"
#include "parameterization-bilinear.hh"
#include "ribbon-compatible.hh"
#include "surface-side-based.hh"

using DomainType = DomainRegular;
using ParamType = ParameterizationBilinear;
using RibbonType = RibbonCompatible;

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
  for (size_t i = 0; i < n_; ++i) {
    Point2D sd2(std::min(std::max(sds[i][0], 0.0), 1.0), gamma(sds[i][1]));
    p += ribbons_[i]->eval(sd2) * blends[i];
  }
  return p;
}

std::shared_ptr<Ribbon>
SurfaceSideBased::newRibbon() const {
  return std::make_shared<RibbonType>();
}
