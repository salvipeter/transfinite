#pragma once

#include "geometry.hh"

namespace Transfinite
{

  class Surface;
  class Domain;

  class Parameterization
  {
  public:
    Parameterization(Surface *surface);
    virtual ~Parameterization();
    void setDomain(Domain *new_domain, size_t i);
    void invalidate();
    virtual Point2D mapToRibbon(Point2D const &uv) const = 0;
    Point2DVector mapToRibbons(Point2DVector const &points) const;
  protected:
    Surface *surface;
    Domain *domain;
    size_t index;
    typedef std::map<Point2D, Point2DVector> ParameterizationCache;
    mutable ParameterizationCache cache;
  };

} // Transfinite
