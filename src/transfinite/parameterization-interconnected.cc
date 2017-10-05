#include <cmath>

#include "domain.hh"
#include "parameterization-interconnected.hh"
#include "utilities.hh"

namespace Transfinite {

ParameterizationInterconnected::~ParameterizationInterconnected() {
}

Point2D
ParameterizationInterconnected::mapToRibbon(size_t i, const Point2D &uv) const {
  double s_prev = ParameterizationBilinear::mapToRibbon(prev(i), uv)[0];
  double s_next = ParameterizationBilinear::mapToRibbon(next(i), uv)[0];
  Point2D sd = ParameterizationBilinear::mapToRibbon(i, uv);
  double Hs = hermite(0, sd[0]);
  sd[1] = (1.0 - s_prev) * Hs + s_next * (1.0 - Hs);
  return sd;
}

Point2DVector
ParameterizationInterconnected::mapToRibbons(const Point2D &uv) const {
  Point2DVector sds = ParameterizationBilinear::mapToRibbons(uv);
  for (size_t i = 0; i < n_; ++i) {
    double Hs = hermite(0, sds[i][0]);
    sds[i][1] = (1.0 - sds[prev(i)][0]) * Hs + sds[next(i)][0] * (1.0 - Hs);
  }
  return sds;
}

} // namespace Transfinite
