#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceGeneralizedCoons : public Surface {
public:
  SurfaceGeneralizedCoons();
  SurfaceGeneralizedCoons(const SurfaceGeneralizedCoons &) = default;
  virtual ~SurfaceGeneralizedCoons();
  SurfaceGeneralizedCoons &operator=(const SurfaceGeneralizedCoons &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
