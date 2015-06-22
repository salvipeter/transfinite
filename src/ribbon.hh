#pragma once

#include "geometry.hh"

namespace Transfinite
{

  class Surface;

  class Ribbon
  {
  public:
    Ribbon(Surface *surface);
    virtual ~Ribbon();
    void setCurve(BSCurve *curve);
    void setFence(Fence *fence);
    void invalidate();
    Point3D evaluate(Point2D const &sd) const;
    Point3DVector evaluate(Point2DVector const &points) const;
  protected:
    Surface *surface;
    BSCurve *curve;
    Fence *fence;
  };

} // Transfinite
