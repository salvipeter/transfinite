#pragma once

#include "ribbon-compatible.hh"

namespace Transfinite {

class RibbonCompatibleWithHandler : public RibbonCompatible {
public:
  virtual ~RibbonCompatibleWithHandler();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const;

protected:
  Vector3D central_;
};

} // namespace Transfinite
