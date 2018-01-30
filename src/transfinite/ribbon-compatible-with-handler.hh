#pragma once

#include "ribbon-compatible.hh"

namespace Transfinite {

class RibbonCompatibleWithHandler : public RibbonCompatible {
public:
  virtual ~RibbonCompatibleWithHandler();
  virtual void update() override;
  virtual Vector3D crossDerivative(double s) const override;

protected:
  Vector3D central_;
};

} // namespace Transfinite
