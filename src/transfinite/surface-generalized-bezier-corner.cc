#include <cassert>

#include "domain-circle.hh"
#include "parameterization-circular.hh"
// #include "domain-regular.hh"
// #include "parameterization-constrained-barycentric.hh"
#include "ribbon-dummy.hh"
#include "surface-generalized-bezier-corner.hh"
#include "utilities.hh"

namespace Transfinite {

using DomainType = DomainCircle;
using ParamType = ParameterizationCircular;
// using DomainType = DomainRegular;
// using ParamType = ParameterizationConstrainedBarycentric;
using RibbonType = RibbonDummy;

SurfaceGeneralizedBezierCorner::SurfaceGeneralizedBezierCorner() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceGeneralizedBezierCorner::~SurfaceGeneralizedBezierCorner() {
}

void
SurfaceGeneralizedBezierCorner::initNetwork(size_t n, size_t degree) {
  assert(degree % 2 == 1 && "This representation works only for odd degrees");
  SurfaceGeneralizedBezier::initNetwork(n, degree);
}

Point3D
SurfaceGeneralizedBezierCorner::eval(const Point2D &uv) const {
  Point3D surface_point(0,0,0);
  Point2DVector sds = param_->mapToRibbons(uv);
  double weight_sum = 0.0;
  for (size_t i = 0; i < n_; ++i) {
    const double &di   = sds[i][1];
    const double &di1  = sds[next(i)][1];
    DoubleVector bl_di, bl_di1;
    bernstein(degree_, di, bl_di);
    bernstein(degree_, di1, bl_di1);
    for (size_t j = 0; j < layers_; ++j) {
      for (size_t k = 0; k < layers_; ++k) {
        double blend = bl_di1[j] * bl_di[k];
        surface_point += nets_[i][degree_-j][k] * blend;
        weight_sum += blend;
      }
    }
  }
  surface_point += central_cp_ * (1.0 - weight_sum);
  return surface_point;
}

double
SurfaceGeneralizedBezierCorner::cornerWeight(size_t i, size_t j, size_t k, const Point2D &uv) const
{
  Point2DVector sds = param_->mapToRibbons(uv);
  const double &di   = sds[i][1];
  const double &di1  = sds[next(i)][1];
  DoubleVector bl_di, bl_di1;
  bernstein(degree_, di, bl_di);
  bernstein(degree_, di1, bl_di1);
  return bl_di1[j] * bl_di[k];
}

double
SurfaceGeneralizedBezierCorner::weight(size_t i, size_t j, size_t k, const Point2D &uv) const {
  if (j < layers_)
    return cornerWeight(prev(i), k, j, uv);
  return cornerWeight(i, degree_ - j, k, uv);
}

std::shared_ptr<Ribbon>
SurfaceGeneralizedBezierCorner::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
