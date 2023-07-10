#pragma once

#include "curves.hh"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Notes:
// - side i consists of vertices i-1 and i
// - in the local coordinate system of side i, these vertices are (0,0) and (1,0)

namespace Transfinite {

using namespace Geometry;

class Domain {
public:
  Domain();
  virtual ~Domain();
  void setSide(size_t i, const std::shared_ptr<Curve> &curve);
  void setSides(const CurveVector &curves);
  virtual bool update();
  size_t size() const;
  virtual size_t meshSize(size_t resolution) const;
  const Point2DVector &parameters(size_t resolution) const;
  virtual TriMesh meshTopology(size_t resolution) const;
  virtual bool onEdge(size_t resolution, size_t index) const;
  const Point2D &center() const;
  virtual Point2D edgePoint(size_t i, double s) const;
  double edgeLength(size_t i) const;
  double angle(size_t i) const;
  const Point2DVector &vertices() const;
  Point2D toLocal(size_t i, const Vector2D &v) const;
  Point2D fromLocal(size_t i, const Vector2D &v) const;
  bool intersectEdgeWithRay(size_t i, const Point2D &p, const Vector2D &v, Point2D &result) const;

protected:
  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }
  virtual void computeCenter();
  virtual Point2DVector parametersImpl(size_t resolution) const;

  CurveVector curves_;
  size_t n_;
  Point2D center_;
  Point2DVector vertices_;
  Vector2DVector du_, dv_;

private:
  mutable Point2DVector parameters_; // cache
};

} // namespace Transfinite
