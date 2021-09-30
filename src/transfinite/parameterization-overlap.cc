#include <cassert>

#include "parameterization-overlap.hh"

namespace Transfinite {

ParameterizationOverlap::~ParameterizationOverlap() {
}

Point2D ParameterizationOverlap::mapToRibbon(size_t i, const Point2D &uv) const {
  assert(n_ % 2 == 0);
  auto bc = barycentric(uv);
  double result = 0.0;
  for (size_t j = 0; j < n_ / 2; ++j)
    result += bc[(i+j)%n_];
  return { result, 0.0 };
}

} // namespace Transfinite
