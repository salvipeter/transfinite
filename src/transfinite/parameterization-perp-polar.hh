#pragma once

#include "parameterization-barycentric.hh"

namespace Transfinite {

class ParameterizationPerpPolar : public Parameterization {
public:
  virtual ~ParameterizationPerpPolar();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};

} // namespace Transfinite
