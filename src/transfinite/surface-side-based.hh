#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceSideBased : public Surface {
public:
  SurfaceSideBased();
  SurfaceSideBased(const SurfaceSideBased &) = default;
  virtual ~SurfaceSideBased();
  SurfaceSideBased &operator=(const SurfaceSideBased &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
