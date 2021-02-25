#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-hybrid.hh"
#include "utilities.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceHybrid::SurfaceHybrid() {
  useSquaredRationalWeights(true);
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceHybrid::~SurfaceHybrid() {
}

Point3D
SurfaceHybrid::eval(const Point2D &uv) const {
  Point3D surface_point(0,0,0);
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends(n_, 0.0);

  double weight_sum = 0.0;
  for (size_t i = 0; i < n_; ++i)
    for (size_t k = 0; k < layers_; ++k)
      for (size_t j = 0; j <= degree_; ++j) {
        double blend = weight(i, j, k, uv);
        if (k >= 2)
          surface_point += nets_[i][j][k] * blend;
        else
          blends[i] += blend;
        weight_sum += blend;
      }
  surface_point += central_cp_ * (1.0 - weight_sum);

  for (size_t i = 0; i < n_; ++i)
    surface_point += sideInterpolant(i, sds[i][0], sds[i][1]) * blends[i];

  return surface_point;
}

std::shared_ptr<Ribbon>
SurfaceHybrid::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
