#pragma once

#include "geometry.hh"

// Notes:
// - side i consists of vertices i-1 and i
// - in the local coordinate system of side i, these vertices are (0,0) and (1,0)

class Domain {
public:
  Domain();
  virtual ~Domain();
  void setSide(size_t i, const std::shared_ptr<BSCurve> &curve);
  void setSides(const CurveVector &curves);
  virtual bool update();
  const Point2DVector &parameters(size_t resolution) const;
  TriMesh meshTopology(size_t resolution) const;
  const Point2D &center() const;
  double edgeLength(size_t i) const;
  double angle(size_t i) const;
  const Point2DVector &vertices() const;
  Point2D toLocal(size_t i, const Vector2D &v) const;

protected:
  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }
  virtual void computeCenter() = 0;

  CurveVector curves_;
  size_t n_;
  Point2D center_;
  Point2DVector vertices_;
  Vector2DVector du_, dv_;

private:
  mutable Point2DVector parameters_; // cache
};
