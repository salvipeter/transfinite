#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceGeneralizedBezier : public Surface {
public:
  SurfaceGeneralizedBezier();
  SurfaceGeneralizedBezier(const SurfaceGeneralizedBezier &) = default;
  virtual ~SurfaceGeneralizedBezier();
  SurfaceGeneralizedBezier &operator=(const SurfaceGeneralizedBezier &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
  size_t degree() const;
  size_t layers() const;
  void initNetwork(size_t n, size_t degree);
  virtual void setupLoop() override;
  Point3D centralControlPoint() const;
  void setCentralControlPoint(const Point3D &p);
  Point3D controlPoint(size_t i, size_t j, size_t k) const;
  void setControlPoint(size_t i, size_t j, size_t k, const Point3D &p);
  double weight(size_t i, size_t j, size_t k, const Point2D &uv) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

private:
  using ControlNet = std::vector<PointVector>;

  size_t degree_, layers_;
  Point3D central_cp_;
  std::vector<ControlNet> nets_;
};

} // namespace Transfinite
