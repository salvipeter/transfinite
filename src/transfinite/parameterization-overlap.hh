#pragma once

#include "parameterization-barycentric.hh"

namespace Transfinite {

class ParameterizationOverlap : public ParameterizationBarycentric {
public:
  virtual ~ParameterizationOverlap();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};

} // namespace Transfinite
