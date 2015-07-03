#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceGeneralizedCoons : public Surface {
public:
  SurfaceGeneralizedCoons();
  virtual ~SurfaceGeneralizedCoons();
  virtual Point3D eval(const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
};

} // namespace Transfinite
