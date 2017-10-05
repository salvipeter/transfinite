#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonNSided : public Ribbon {
public:
  virtual ~RibbonNSided();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const;
  virtual Point3D eval(const Point2D &sd) const;
  virtual Vector3D normal(double s) const;

protected:
  double base_length_;
};

} // namespace Transfinite
