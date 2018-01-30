#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfacePolar : public Surface {
public:
  SurfacePolar(size_t max_degree = 5);
  SurfacePolar(const SurfacePolar &) = default;
  virtual ~SurfacePolar();
  SurfacePolar &operator=(const SurfacePolar &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
  Point3D polarRibbon(size_t i, const Point2D &pd) const;

private:
  size_t max_degree_;
};

} // namespace Transfinite
