#pragma once

#include "parameterization.hh"

class ParameterizationBilinear : public Parameterization {
public:
  ParameterizationBilinear(Surface *surface);
  virtual ~ParameterizationBilinear();
  virtual Point2D mapToRibbon(const Point2D &uv) const;
};
