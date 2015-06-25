#pragma once

#include "geometry.hh"

class Surface;

class Ribbon {
public:
  Ribbon(Surface *surface);
  virtual ~Ribbon();
  void setCurve(BSCurve *curve);
  void setFence(Fence *fence);
  void invalidate();
  Point3D evaluate(const Point2D &sd) const;
  Point3DVector evaluate(const Point2DVector &points) const;
protected:
  Surface *surface_;
  BSCurve *curve_;
  Fence *fence_;
};
