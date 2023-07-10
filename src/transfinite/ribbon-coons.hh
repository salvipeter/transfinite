#pragma once

#include "ribbon.hh"

namespace Transfinite {

class RibbonCoons : public Ribbon {
public:
  virtual ~RibbonCoons();
  virtual void update() override;
  virtual Vector3D crossDerivative(double s) const override;
  virtual Point3D eval(const Point2D &sd) const override;

protected:
  std::shared_ptr<Curve> left_, right_, top_;
  Point3D bl_, br_, tl_, tr_;
};

} // namespace Transfinite
