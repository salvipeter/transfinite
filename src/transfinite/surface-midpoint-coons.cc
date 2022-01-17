#include <algorithm>
#include <numeric>

#include "domain-regular.hh"
#include "parameterization-constrained-barycentric.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-midpoint-coons.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationConstrainedBarycentric;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceMidpointCoons::SurfaceMidpointCoons() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceMidpointCoons::~SurfaceMidpointCoons() {
}

Point3D
SurfaceMidpointCoons::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCornerDeficient(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i) {
    double s = sds[i][0], d = sds[i][1], s1 = sds[next(i)][0];
    p += sideInterpolant(i, s, d) * (blends[i] + blends[prev(i)]);
    p -= cornerCorrection(i, 1.0 - s, s1) * blends[i];
  }
  p += central_cp_ * (1.0 - std::accumulate(blends.begin(), blends.end(), 0.0));
  return p;
}

} // namespace Transfinite
