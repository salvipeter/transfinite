#include "domain-circular.hh"

#include <algorithm>
#include <numeric>

namespace Transfinite {

DomainCircular::~DomainCircular() {
}

bool
DomainCircular::update() {
  n_ = curves_.size();

  DoubleVector lengths; lengths.reserve(n_);
  std::transform(curves_.begin(), curves_.end(), std::back_inserter(lengths),
                 [](const std::shared_ptr<Curve> &c) { return c->arcLength(0.0, 1.0); });
  double normalizer = 2.0 * M_PI / std::accumulate(lengths.begin(), lengths.end(), 0.0);
  std::transform(lengths.begin(), lengths.end(), lengths.begin(),
                 [normalizer](double x) { return x * normalizer; });

  double alpha = 0.0;
  vertices_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    alpha += lengths[i];
    vertices_[i] = Point2D(std::cos(alpha), std::sin(alpha));
  }

  return Domain::update();
}

} // namespace Transfinite
