#include <cmath>

#include "domain.hh"
#include "parameterization-parallel.hh"

namespace Transfinite {

ParameterizationParallel::~ParameterizationParallel() {
}

Point2D
ParameterizationParallel::mapToRibbon(size_t i, const Point2D &uv) const {
  Point2D sd = ParameterizationBilinear::mapToRibbon(i, uv);
  sd[1] = domain_->toLocal(i, uv - domain_->vertices()[i])[1] * multipliers_[i];
  return sd;
}

void
ParameterizationParallel::updateMultipliers() {
  multipliers_.resize(n_);
  const Point2DVector &v = domain_->vertices();
  for (size_t i = 0; i < n_; ++i) {
    double max_dist = 0.0;
    for (const auto &p : v) {
      double dist = domain_->toLocal(i, p - v[i])[1];
      if (dist > max_dist)
        max_dist = dist;
    }
    multipliers_[i] = 1.0 / max_dist;
  }
}

void
ParameterizationParallel::update() {
  // Should be called first to update n_
  ParameterizationBilinear::update();
  updateMultipliers();
}

} // namespace Transfinite
