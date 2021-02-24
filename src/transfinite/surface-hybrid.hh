#pragma once

#include "surface-generalized-bezier.hh"

namespace Transfinite {

class SurfaceHybrid : public SurfaceGeneralizedBezier {
public:
  SurfaceHybrid();
  SurfaceHybrid(const SurfaceHybrid &) = default;
  virtual ~SurfaceHybrid();
  SurfaceHybrid &operator=(const SurfaceHybrid &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
