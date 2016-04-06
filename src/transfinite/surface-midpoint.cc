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
  // Reinitialize parameterization, because the parent class used ParameterizationBilinear
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceMidpoint::~SurfaceMidpoint() {
}

void
SurfaceMidpoint::update(size_t i) {
  SurfaceCornerBased::update(i);
  updateCentralControlPoint();
}

void
SurfaceMidpoint::update() {
  SurfaceCornerBased::update();
  updateCentralControlPoint();
}

Point3D
SurfaceMidpoint::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = deficientCornerBlend(sds);
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
SurfaceMidpoint::hermite(double u) {
  double u1 = 1.0 - u, u1_2 = u1 * u1;
  return u1_2 * u1 + 3.0 * u1_2 * u;
}

DoubleVector
SurfaceMidpoint::deficientCornerBlend(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    size_t ip = next(i);
    if (sds[i][1] < epsilon && sds[ip][1] < epsilon) {
      blf.push_back(1.0);
      continue;
    }
    blf.push_back((sds[ip][1] * hermite(1.0 - sds[i][0]) * hermite(sds[i][1] ) +
                   sds[i][1]  * hermite(   sds[ip][0]  ) * hermite(sds[ip][1])) /
                  (sds[i][1] + sds[ip][1]));
  }
  return blf;
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
  DoubleVector blends = deficientCornerBlend(param_->mapToRibbons(center));
  double blf_sum = std::accumulate(blends.begin(), blends.end(), 0.0);
  if (fabs(1.0 - blf_sum) < epsilon)
    central_cp_ = midpoint_;  // as good as anything else
  else
    central_cp_ = (midpoint_ - s) / (1.0 - blf_sum);
}

} // namespace Transfinite
