#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceC0Coons : public Surface {
public:
  SurfaceC0Coons();
  SurfaceC0Coons(const SurfaceC0Coons &) = default;
  virtual ~SurfaceC0Coons();
  SurfaceC0Coons &operator=(const SurfaceC0Coons &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
