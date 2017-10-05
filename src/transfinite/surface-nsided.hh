#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceNSided : public Surface {
public:
  SurfaceNSided();
  SurfaceNSided(const SurfaceNSided &) = default;
  virtual ~SurfaceNSided();
  SurfaceNSided &operator=(const SurfaceNSided &) = default;
  virtual void update(size_t i);
  virtual void update();
  virtual Point3D eval(const Point2D &uv) const;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;

  std::shared_ptr<Parameterization> blend_param_;
};

} // namespace Transfinite
