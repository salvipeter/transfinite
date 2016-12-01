#pragma once

#include "surface-corner-based.hh"

namespace Transfinite {

class SurfaceMidpoint : public SurfaceCornerBased {
public:
  SurfaceMidpoint();
  SurfaceMidpoint(const SurfaceMidpoint &) = default;
  virtual ~SurfaceMidpoint();
  SurfaceMidpoint &operator=(const SurfaceMidpoint &) = default;
  virtual void update(size_t i);
  virtual void update();
  virtual Point3D eval(const Point2D &uv) const;
  using Surface::eval;
  void setMidpoint(const Point3D &p);
  void unsetMidpoint();

private:
  void updateCentralControlPoint();

  bool midpoint_set_;
  Point3D midpoint_, central_cp_;
};

} // namespace Transfinite
