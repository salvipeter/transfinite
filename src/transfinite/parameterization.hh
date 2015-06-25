#pragma once

#include "geometry.hh"

class Surface;
class Domain;

class Parameterization {
public:
  Parameterization(Surface *surface);
  virtual ~Parameterization();
  void setDomain(Domain *new_domain, size_t i);
  void invalidate();
  virtual Point2D mapToRibbon(const Point2D &uv) const = 0;
  Point2DVector mapToRibbons(const Point2DVector &points) const;

protected:
  Surface *surface_;
  Domain *domain_;
  size_t index_;
  typedef std::map<Point2D, Point2DVector> ParameterizationCache;
  mutable ParameterizationCache cache_;
};
