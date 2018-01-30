#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonDummy : public Ribbon {
public:
  virtual ~RibbonDummy() {}
  virtual void update() override {}
  virtual Vector3D crossDerivative(double) const override { return Vector3D(0,0,0); }
};

} // namespace Transfinite
