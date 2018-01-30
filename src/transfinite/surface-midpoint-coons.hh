#pragma once

#include "surface-midpoint.hh"

namespace Transfinite {

class SurfaceMidpointCoons : public SurfaceMidpoint {
public:
  SurfaceMidpointCoons();
  SurfaceMidpointCoons(const SurfaceMidpointCoons &) = default;
  virtual ~SurfaceMidpointCoons();
  SurfaceMidpointCoons &operator=(const SurfaceMidpointCoons &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
};

} // namespace Transfinite
