#pragma once

#include "surface.hh"

class SurfaceCornerBased : public Surface {
public:
  SurfaceCornerBased();
  virtual ~SurfaceCornerBased();
  virtual Point3D eval(const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
  Point3D cornerInterpolant(size_t i, const Point2DVector &sds) const;
};
