#include <cmath>

#include "domain.hh"
#include "parameterization-superd.hh"

namespace Transfinite {

ParameterizationSuperD::~ParameterizationSuperD() {
}

void
ParameterizationSuperD::updateMultipliers() { // assumes regular domain
  const Point2DVector &v = domain_->vertices();
  // In SuperD, currently k = 1 for n <= 8 and n/2 for n > 8
  size_t k = n_ <= 6 ? 1 : (n_ - 3) / 2; // same as std::max(1, (n_ - 3) / 2)
  double max_distance = domain_->toLocal(0, (v[k] - v[0]))[1];
  multipliers_.assign(n_, 1.0 / max_distance);
}

} // namespace Transfinite
