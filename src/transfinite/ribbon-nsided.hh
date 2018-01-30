#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonNSided : public Ribbon {
public:
  virtual ~RibbonNSided();
  virtual void update() override;
  virtual Vector3D crossDerivative(double s) const override;
  virtual Point3D eval(const Point2D &sd) const override;
  virtual Vector3D normal(double s) const override;

protected:
  double base_length_;
};

} // namespace Transfinite
