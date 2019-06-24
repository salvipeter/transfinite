#pragma once

#include "parameterization.hh"

namespace Transfinite {

class ParameterizationBarycentric : public Parameterization {
public:
  enum class BarycentricType { WACHSPRESS, MEAN_VALUE, HARMONIC };

  ParameterizationBarycentric();
  ParameterizationBarycentric(BarycentricType type);
  virtual ~ParameterizationBarycentric();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
  virtual void update() override;
  const DoubleVector &barycentric(const Point2D &uv) const;

private:
  const BarycentricType type_;

  using CoordinateCache = std::map<Point2D, DoubleVector, PointComparator>;
  mutable CoordinateCache cache_;
};

} // namespace Transfinite
