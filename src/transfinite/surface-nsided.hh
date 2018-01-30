#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceNSided : public Surface {
public:
  SurfaceNSided();
  SurfaceNSided(const SurfaceNSided &) = default;
  virtual ~SurfaceNSided();
  SurfaceNSided &operator=(const SurfaceNSided &) = default;
  virtual void update(size_t i) override;
  virtual void update() override;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

  std::shared_ptr<Parameterization> blend_param_;
};

} // namespace Transfinite
