#pragma once

#include "parameterization.hh"

namespace Transfinite {

// h-only parameterization (returns same value in both coordinates)

class ParameterizationCircular : public Parameterization {
public:
  virtual ~ParameterizationCircular();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};

}
