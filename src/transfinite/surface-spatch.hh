#pragma once

#include <map>

#include "surface.hh"

namespace Transfinite {

class SurfaceSPatch : public Surface {
public:
  using Index = std::vector<size_t>;
  SurfaceSPatch();
  SurfaceSPatch(const SurfaceSPatch &) = default;
  virtual ~SurfaceSPatch();
  SurfaceSPatch &operator=(const SurfaceSPatch &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
  void initNetwork(size_t n, size_t d);
  virtual void setupLoop() override;
  void setControlPoint(const Index &i, const Point3D &p);
  Point3D controlPoint(const Index &i) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

private:
  size_t depth_;
  std::map<Index, Point3D> net_;
};

} // namespace Transfinite
