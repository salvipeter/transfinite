#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceCompositeRibbon : public Surface {
public:
  SurfaceCompositeRibbon();
  virtual ~SurfaceCompositeRibbon();
  virtual Point3D eval(const Point2D &uv) const;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
  Point3D compositeRibbon(size_t i, const Point2DVector &sds) const;
};

} // namespace Transfinite
