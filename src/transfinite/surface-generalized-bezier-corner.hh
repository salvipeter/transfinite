#pragma once

#include "surface-generalized-bezier.hh"

namespace Transfinite {

class SurfaceGeneralizedBezierCorner : public SurfaceGeneralizedBezier {
public:
  SurfaceGeneralizedBezierCorner();
  SurfaceGeneralizedBezierCorner(const SurfaceGeneralizedBezierCorner &) = default;
  virtual ~SurfaceGeneralizedBezierCorner();
  virtual void initNetwork(size_t n, size_t degree) override;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
  virtual double weight(size_t i, size_t j, size_t k, const Point2D &uv) const override;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

private:
  double cornerWeight(size_t i, size_t j, size_t k, const Point2D &uv) const;
};

} // namespace Transfinite
