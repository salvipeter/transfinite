#include "domain-circular.hh"

#include <algorithm>
#include <cmath>

DomainCircular::~DomainCircular() {
}

void
DomainCircular::setSides(const CurveVector &curves) {
  size_t m = curves.size();

  DoubleVector lengths;
  lengths.reserve(curves.size());
  std::transform(curves.begin(), curves.end(), std::back_inserter(lengths),
                 [](const BSCurve &c) { return c.arcLength(0.0, 1.0); });
  double normalizer = 2.0 * M_PI / std::accumulate(lengths.begin(), lengths.end(), 0.0);
  std::transform(lengths.begin(), lengths.end(), lengths.begin(),
                 [normalizer](double x) { return x * normalizer; });

  double alpha = 0.0;
  vertices_.resize(m);
  for(size_t i = 0; i < m; ++i) {
    alpha += lengths[i];
    vertices_[i] = Point2D(std::cos(alpha), std::sin(alpha));
  }

  invalidate();
}

void
DomainCircular::computeCenter() {
  DoubleVector lengths;
  lengths.reserve(n_);
  for (size_t i = 0; i < n_; ++i)
    lengths.push_back((vertices_[i] - vertices_[prev(i)]).norm());
  center_ = Point2D(0, 0);
  for (size_t i = 0; i < n_; ++i)
    center_ += vertices_[i] * (lengths[i] + lengths[next(i)]);
  center_ /= std::accumulate(lengths.begin(), lengths.end(), 0.0) * 2.0;
}
