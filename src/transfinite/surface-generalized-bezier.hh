#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceGeneralizedBezier : public Surface {
public:
  SurfaceGeneralizedBezier();
  virtual ~SurfaceGeneralizedBezier();
  virtual Point3D eval(const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
};

} // namespace Transfinite
