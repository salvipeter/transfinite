#include "parameterization.hh"

Parameterization::Parameterization(Surface *surface)
  : surface(surface) {
}

Parameterization::~Parameterization() {
}

void
Parameterization::setDomain(Domain *new_domain, size_t i) {
  domain = new_domain;
  index = i;
  invalidate();
}

void
Parameterization::invalidate() {
  cache.clear();
}

Point2DVector
Parameterization::mapToRibbon(const Point2DVector &points) const {
  Point2DVector result;
  result.reserve(points.size());
  for(Point2DVector::const_iterator i = points.begin(), ie = points.end(); i != ie; ++i)
    result.push_back(mapToRibbon(*i));
  return result;
}
