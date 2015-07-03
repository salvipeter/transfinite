#pragma once

#include "parameterization-bilinear.hh"

namespace Transfinite {

class ParameterizationInterconnected : public ParameterizationBilinear {
public:
  virtual ~ParameterizationInterconnected();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const;
  virtual Point2DVector mapToRibbons(const Point2D &uv) const;

protected:
  static double blendHermite(double x);
};

} // namespace Transfinite
