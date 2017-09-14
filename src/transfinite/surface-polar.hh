#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfacePolar : public Surface {
public:
  SurfacePolar();
  SurfacePolar(const SurfacePolar &) = default;
  virtual ~SurfacePolar();
  SurfacePolar &operator=(const SurfacePolar &) = default;
  virtual Point3D eval(const Point2D &uv) const;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
  Point3D polarRibbon(size_t i, const Point2DVector &pds) const;
};

} // namespace Transfinite
