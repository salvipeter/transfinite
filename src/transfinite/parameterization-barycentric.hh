#pragma once

#include "parameterization.hh"

class ParameterizationBarycentric : public Parameterization {
public:
  virtual ~ParameterizationBarycentric();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const;
  virtual void invalidate();

protected:
  const DoubleVector &barycentric(const Point2D &uv) const;

private:
  enum class BarycentricType { WACHSPRESS, MEAN_VALUE, HARMONIC };
  static const BarycentricType type_ = BarycentricType::WACHSPRESS;

  using CoordinateCache = std::map<Point2D, DoubleVector, PointComparator>;
  mutable CoordinateCache cache_;
};
