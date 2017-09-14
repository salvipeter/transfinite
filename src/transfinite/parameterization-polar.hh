#pragma once

#include "parameterization-barycentric.hh"

namespace Transfinite {

class ParameterizationPolar : public ParameterizationBarycentric {
public:
  ParameterizationPolar();
  virtual ~ParameterizationPolar();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const;
};

} // namespace Transfinite
