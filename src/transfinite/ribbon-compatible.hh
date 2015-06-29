#pragma once

#include "ribbon.hh"

class RibbonCompatible : public Ribbon {
public:
  virtual ~RibbonCompatible();
  virtual void update();
  virtual Point3D eval(const Point2D &sd) const;

protected:
  Vector3D prev_tangent_, next_tangent_;
};
