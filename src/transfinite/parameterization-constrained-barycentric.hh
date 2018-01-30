#pragma once

#include "parameterization-barycentric.hh"

namespace Transfinite {

class ParameterizationConstrainedBarycentric : public ParameterizationBarycentric {
public:
  virtual ~ParameterizationConstrainedBarycentric();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};

} // namespace Transfinite
