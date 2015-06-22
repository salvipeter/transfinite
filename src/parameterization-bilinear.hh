#pragma once

#include "parameterization.hh"

namespace Transfinite
{

  class ParameterizationBilinear : public Parameterization
  {
  public:
    ParameterizationBilinear(Surface *surface);
    virtual ~ParameterizationBilinear();
    virtual Point2D mapToRibbon(Point2D const &uv) const;
  };

} // Transfinite
