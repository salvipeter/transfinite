#include <algorithm>
#include <stdexcept>

#include "domain.hh"
#include "parameterization.hh"

namespace Transfinite {

Parameterization::~Parameterization() {
}

void
Parameterization::setDomain(const std::shared_ptr<Domain> &new_domain) {
  domain_ = new_domain;
}

void
Parameterization::update() {
  n_ = domain_->vertices().size();
  cache_.clear();
}

Point2DVector
Parameterization::mapToRibbons(const Point2D &uv) const {
  auto cached = cache_.find(uv);
  if (cached != cache_.end())
    return cached->second;

  Point2DVector result; result.reserve(n_);
  for (size_t i = 0; i < n_; ++i)
    result.push_back(mapToRibbon(i, uv));

  cache_[uv] = result;
  return result;
}

Point2D
Parameterization::inverse(size_t i, const Point2D &pd) const {
  throw std::logic_error("inverse() is not implemented for this parameterization");
}

} // namespace Transfinite
