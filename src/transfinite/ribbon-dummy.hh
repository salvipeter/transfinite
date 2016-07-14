#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonDummy : public Ribbon {
public:
  virtual ~RibbonDummy() {}
  virtual void update() {}
  virtual Vector3D crossDerivative(double) const { return Vector3D(0,0,0); }
};

} // namespace Transfinite
