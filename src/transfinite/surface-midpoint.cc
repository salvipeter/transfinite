#include <algorithm>
#include <numeric>

#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-midpoint.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceMidpoint::SurfaceMidpoint() : midpoint_set_(false) {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceMidpoint::~SurfaceMidpoint() {
}

void
SurfaceMidpoint::update(size_t i) {
  Surface::update(i);
  updateCentralControlPoint();
}

void
SurfaceMidpoint::update() {
  Surface::update();
  updateCentralControlPoint();
}

Point3D
SurfaceMidpoint::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCornerDeficient(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += cornerInterpolant(i, sds) * blends[i];
  p += central_cp_ * (1.0 - std::accumulate(blends.begin(), blends.end(), 0.0));
  return p;
}

void
SurfaceMidpoint::setMidpoint(const Point3D &p) {
  midpoint_ = p;
  midpoint_set_ = true;
  updateCentralControlPoint();
}
 
void
SurfaceMidpoint::unsetMidpoint() {
  midpoint_set_ = false;
  updateCentralControlPoint();
}

double
SurfaceMidpoint::deficiency(const Point2D &p) const {
  DoubleVector blends = blendCornerDeficient(param_->mapToRibbons(p));
  double blf_sum = std::accumulate(blends.begin(), blends.end(), 0.0);
  return 1.0 - blf_sum;
}

std::shared_ptr<Ribbon>
SurfaceMidpoint::newRibbon() const {
  return std::make_shared<RibbonType>();
}

void
SurfaceMidpoint::updateCentralControlPoint() {
  if (!midpoint_set_) {
    // Compute a default midpoint position
    midpoint_ = Point3D(0,0,0);
    for (size_t i = 0; i < n_; ++i)
      midpoint_ += sideInterpolant(i, 0.5, 0.5);
    midpoint_ /= n_;
  }
  Point2D center = domain_->center();
  central_cp_ = Point3D(0,0,0);
  Point3D s = eval(center);
  double def = deficiency(center);
  if (fabs(def) < epsilon)
    central_cp_ = midpoint_;  // as good as anything else
  else
    central_cp_ = (midpoint_ - s) / def;
}

} // namespace Transfinite
