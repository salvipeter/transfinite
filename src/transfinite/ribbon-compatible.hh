#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonCompatible : public Ribbon {
public:
  virtual ~RibbonCompatible();
  virtual void update() override;
  virtual Vector3D crossDerivative(double s) const override;

protected:
  Vector3D prev_tangent_, next_tangent_;
};

} // namespace Transfinite
