#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceBiharmonic : public Surface {
public:
  SurfaceBiharmonic();
  SurfaceBiharmonic(const SurfaceBiharmonic &) = default;
  virtual ~SurfaceBiharmonic();
  SurfaceBiharmonic &operator=(const SurfaceBiharmonic &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  virtual TriMesh eval(size_t resolution) const override;

protected:
  auto generateDomainOld(size_t resolution) const;
  auto generateDomain(size_t resolution) const;
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};

} // namespace Transfinite
