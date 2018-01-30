#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceCompositeRibbon : public Surface {
public:
  SurfaceCompositeRibbon();
  SurfaceCompositeRibbon(const SurfaceCompositeRibbon &) = default;
  virtual ~SurfaceCompositeRibbon();
  SurfaceCompositeRibbon &operator=(const SurfaceCompositeRibbon &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
  Point3D compositeRibbon(size_t i, const Point2D &sd) const;
};

} // namespace Transfinite
