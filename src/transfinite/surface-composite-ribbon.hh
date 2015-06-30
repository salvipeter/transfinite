#pragma once

#include "surface.hh"

class SurfaceCompositeRibbon : public Surface {
public:
  SurfaceCompositeRibbon();
  virtual ~SurfaceCompositeRibbon();
  virtual Point3D eval(const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
  Point3D compositeRibbon(size_t i, const Point2DVector &sds) const;
};
