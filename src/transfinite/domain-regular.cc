#include "domain-regular.hh"

#include <cmath>

DomainRegular::~DomainRegular() {
}

void
DomainRegular::setSides(const CurveVector &curves) {
  size_t m = curves.size();
  if(n_ == m)
    return;

  double alpha = 2.0 * M_PI / m;
  vertices_.resize(m);
  for(size_t i = 0; i < m; ++i)
    vertices_[i] = Point2D(std::cos(alpha * i), std::sin(alpha * i));

  invalidate();
}

void
DomainRegular::computeCenter() {
  center_ = Point2D(0.0, 0.0);
}
