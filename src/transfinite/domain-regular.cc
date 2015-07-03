#include "domain-regular.hh"

#include <cmath>

namespace Transfinite {

DomainRegular::~DomainRegular() {
}

bool
DomainRegular::update() {
  size_t m = curves_.size();
  if(n_ == m)
    return false;

  double alpha = 2.0 * M_PI / m;
  vertices_.resize(m);
  for(size_t i = 0; i < m; ++i)
    vertices_[i] = Point2D(std::cos(alpha * i), std::sin(alpha * i));

  return Domain::update();
}

void
DomainRegular::computeCenter() {
  center_ = Point2D(0.0, 0.0);
}

} // namespace Transfinite
