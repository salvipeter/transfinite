#pragma once

#include "geometry.hh"

namespace Transfinite {

class WithHandler {
public:
  WithHandler();
  virtual ~WithHandler();
  void setMultiplier(double m);
  void setHandler(const Vector3D &h);
  virtual void resetHandler();

protected:
  Vector3D handler_;
  double multiplier_;
  bool handler_initialized_;
};

} // namespace Transfinite
