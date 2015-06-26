#include <algorithm>

#include "domain.hh"
#include "parameterization-barycentric.hh"

ParameterizationBarycentric::~ParameterizationBarycentric() {
}

Point2D
ParameterizationBarycentric::mapToRibbon(size_t i, const Point2D &uv) const {
  const DoubleVector &l = barycentric(uv);
  Point2D sd;
  double denom = l[prev(i)] + l[i];
  if (denom < epsilon) {
    sd[0] = 0.5;                // kutykurutty
  } else
    sd[0] = l[i] / denom;
  sd[1] = 1.0 - l[i] - l[prev(i)];
  return sd;
}

void
ParameterizationBarycentric::invalidate() {
  cache_.clear();
  Parameterization::invalidate();
}

const DoubleVector &
ParameterizationBarycentric::barycentric(const Point2D &uv) const {
  auto cached = cache_.find(uv);
  if(cached != cache_.end())
    return cached->second;

  Vector2DVector vectors; vectors.reserve(n_);
  std::transform(domain_->vertices().begin(), domain_->vertices().end(),
                 std::back_inserter(vectors),
                 [uv](const Point2D &p) { return uv - p; });

  DoubleVector areas; areas.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    const Vector2D &si = vectors[i];
    const Vector2D &si1 = vectors[next(i)];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n_);

  for (size_t i = 0; i < n_; ++i) {
    size_t i_1 = prev(i), i1 = next(i);
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n_; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector2D &si_1 = vectors[i_1];
    const Vector2D &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    double ri_1, ri, ri1;
    switch (type_) {
    case BarycentricType::WACHSPRESS:
      ri_1 = 1.0;
      ri = 1.0;
      ri1 = 1.0;
      break;
    case BarycentricType::MEAN_VALUE:
      ri_1 = vectors[i_1].norm();
      ri = vectors[i].norm();
      ri1 = vectors[i1].norm();
      break;
    case BarycentricType::HARMONIC:
      ri_1 = vectors[i_1].normSqr();
      ri = vectors[i].normSqr();
      ri1 = vectors[i1].normSqr();
      break;
    };
    l.push_back(ri_1 * Ai_1 + ri1 * Ai - ri * Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });

  cache_[uv] = l;
  return cache_[uv];
}
