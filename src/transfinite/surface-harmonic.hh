#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceHarmonic : public Surface {
public:
  SurfaceHarmonic();
  SurfaceHarmonic(const SurfaceHarmonic &) = default;
  virtual ~SurfaceHarmonic();
  SurfaceHarmonic &operator=(const SurfaceHarmonic &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  virtual TriMesh eval(size_t resolution) const override;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
