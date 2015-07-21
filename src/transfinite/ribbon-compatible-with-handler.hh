#pragma once

#include "ribbon-compatible.hh"
#include "with-handler.hh"

namespace Transfinite {

class RibbonCompatibleWithHandler : public RibbonCompatible, public WithHandler {
public:
  virtual ~RibbonCompatibleWithHandler();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const;
  virtual void resetHandler();

protected:
  Vector3D central_;
};

} // namespace Transfinite
