#pragma once

#include "surface.hh"

class SurfaceSideBased : public Surface {
public:
  SurfaceSideBased();
  virtual ~SurfaceSideBased();
  virtual Point3D eval(const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const;
};
