#include "domain-angular.hh"
#include "parameterization-parallel.hh"
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
  blend_param_ = std::make_shared<ParameterizationParallel>();
  blend_param_->setDomain(domain_);
}

SurfaceNSided::~SurfaceNSided() {
}

void
SurfaceNSided::update(size_t i) {
  Surface::update(i);
  blend_param_->update();
}

void
SurfaceNSided::update() {
  Surface::update();
  blend_param_->update();
}

Point3D
SurfaceNSided::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  Point2DVector blend_sds = blend_param_->mapToRibbons(uv);
  DoubleVector blends = blendSideSingular(blend_sds);
  Point3D p(0,0,0);
  return ribbons_[3]->eval(sds[3]);
  for (size_t i = 0; i < n_; ++i)
    p += ribbons_[i]->eval(sds[i]) * blends[i];
  return p;
}

std::shared_ptr<Ribbon>
SurfaceNSided::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
