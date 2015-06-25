#pragma once

#include "geometry.hh"

// Notes:
// - side i consists of vertices i-1 and i
// - in the local coordinate system of side i, these vertices are (0,0) and (1,0)

class Domain {
public:
  Domain();
  virtual ~Domain();
  virtual void setSides(const CurveVector &curves) = 0;
  void invalidate();
  const Point2DVector &globalParameters(size_t resolution) const;
  TriMesh meshTopology(size_t resolution) const;
  const Point2D &center() const;
  const Point2DVector &verticesGlobal() const;
  const Point2DVector &verticesLocal(size_t i) const; // in side i's local coordinate system
  Point2D toLocal(size_t i, const Point2D &p) const;
  Point2D toGlobal(size_t i, const Point2D &p) const;

protected:
  size_t next(size_t i) const { return (i + 1) % n_; }
  size_t prev(size_t i) const { return (i + n_ - 1) % n_; }
  virtual void computeCenter() = 0;

  size_t n_;
  Point2D center_;
  Point2DVector vertices_;
  Vector2DVector du_, dv_;
  std::vector<Point2DVector> local_vertices_;
  mutable Point2DVector parameters_; // cache
};
