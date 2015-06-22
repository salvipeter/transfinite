#pragma once

#include "geometry.hh"

// Notes:
// - side i consists of vertices i-1 and i
// - in the local coordinate system of side i, these vertices are (0,0) and (1,0)

namespace Transfinite
{

  class Surface;

  class Domain
  {
  public:
    Domain(Surface *surface);
    virtual ~Domain();
    virtual void setSides(CurveVector const &curves) = 0;
    void invalidate();
    Point2DVector globalParameters(size_t resolution) const;
    Mesh meshTopology(size_t resolution) const;
    virtual Point2D center() const = 0;
    Point2DVector const &verticesGlobal() const;
    Point2DVector const &verticesLocal(size_t i) const; // in side i's local coordinate system
    Point2D toLocal(size_t i, Point2D const &p) const;
    Point2D toGlobal(size_t i, Point2D const &p) const;
  protected:
    size_t next(size_t i) const { return (i + 1) % n; }
    size_t prev(size_t i) const { return (i + n - 1) % n; }

    Surface *surface;
    size_t n;
    Point2DVector vertices;
    std::vector<Vector2D> du, dv;
    std::vector<Point2DVector> local_vertices;
    mutable Point2DVector parameters; // cache
  };

} // Transfinite
