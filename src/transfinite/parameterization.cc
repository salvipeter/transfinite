#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"

Parameterization::~Parameterization() {
}

void
Parameterization::setDomain(const std::shared_ptr<Domain> &new_domain) {
  domain_ = new_domain;
  invalidate();
}

void
Parameterization::invalidate() {
  n_ = domain_->vertices().size();
  cache_.clear();
}

Point2DVector
Parameterization::mapToRibbons(const Point2D &uv) const {
  auto cached = cache_.find(uv);
  if(cached != cache_.end())
    return cached->second;

  Point2DVector result; result.reserve(n_);
  for (size_t i = 0; i < n_; ++i)
    result.push_back(mapToRibbon(i, uv));

  cache_[uv] = result;
  return result;
}
