#pragma once

#include "ribbon.hh"

class RibbonCompatible : public Ribbon {
public:
  virtual ~RibbonCompatible();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const;

protected:
  Vector3D prev_tangent_, next_tangent_;
};
