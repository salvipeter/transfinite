#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonPerpendicular : public Ribbon {
public:
  virtual ~RibbonPerpendicular();
  virtual void update() override;
  virtual Vector3D crossDerivative(double s) const override;

protected:
  double prev_alpha_, prev_beta_, prev_norm_;
  double next_alpha_, next_beta_, next_norm_;
};

} // namespace Transfinite
