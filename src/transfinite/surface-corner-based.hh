#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceCornerBased : public Surface {
public:
  SurfaceCornerBased();
  SurfaceCornerBased(const SurfaceCornerBased &) = default;
  virtual ~SurfaceCornerBased();
  SurfaceCornerBased &operator=(const SurfaceCornerBased &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
