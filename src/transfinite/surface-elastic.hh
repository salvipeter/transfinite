#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceElastic : public Surface {
public:
  SurfaceElastic();
  SurfaceElastic(const SurfaceElastic &) = default;
  virtual ~SurfaceElastic();
  SurfaceElastic &operator=(const SurfaceElastic &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

private:
  double domain_tolerance;
};

} // namespace Transfinite
