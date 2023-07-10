#include "domain-angular.hh"
#include "utilities.hh"

#include <algorithm>
#include <numeric>

namespace Transfinite {

DomainAngular::~DomainAngular() {
}

bool
DomainAngular::update() {
  n_ = curves_.size();

  // Compute lengths
  DoubleVector lengths; lengths.reserve(curves_.size());
  std::transform(curves_.begin(), curves_.end(), std::back_inserter(lengths),
                 [](const std::shared_ptr<Curve> &c) { return c->arcLength(0.0, 1.0); });
  double length_sum = std::accumulate(lengths.begin(), lengths.end(), 0.0);

  // Compute angles
  VectorVector der;
  DoubleVector angles; angles.reserve(n_);
  double angle_sum = 0.0;
  for (size_t i = 0; i < n_; ++i) {
    curves_[i]->eval(1.0, 1, der);
    Vector3D v1 = -der[1].normalize();
    curves_[next(i)]->eval(0.0, 1, der);
    Vector3D v2 = der[1].normalize();
    angles.push_back(std::acos(inrange(-1, v1 * v2, 1)));
    angle_sum += angles.back();
  }
  double angle_multiplier = (n_ - 2) * M_PI / angle_sum;
  std::transform(angles.begin(), angles.end(), angles.begin(),
                 [angle_multiplier](double x) { return M_PI - x * angle_multiplier; });

  // Compute open domain
  vertices_.resize(n_);
  double dir = 0.0;
  Vector2D prev_v(0.0, 0.0);
  for (size_t i = 0; i < n_; ++i) {
    vertices_[i] = prev_v + Vector2D(std::cos(dir), std::sin(dir)) * lengths[i];
    dir += angles[i];
    prev_v = vertices_[i];
  }

  // Compute closed domain
  double accumulated = 0.0;
  for (size_t i = 0; i < n_; ++i) {
    accumulated += lengths[i];
    vertices_[i] -= vertices_.back() * accumulated / length_sum;
  }

  // Rescale to [-1,-1]x[1,1]
  double minx = 0.0, miny = 0.0, maxx = 0.0, maxy = 0.0;
  for (const auto &v : vertices_) {
    minx = std::min(minx, v[0]); miny = std::min(miny, v[1]);
    maxx = std::max(maxx, v[0]); maxy = std::max(maxy, v[1]);
  }
  double width = std::max(maxx - minx, maxy - miny);
  Point2D topleft(-1.0, -1.0), minp(minx, miny);
  for (auto &v : vertices_)
    v = topleft + (v - minp) * 2.0 / width;;

  // Convexity check
  bool convex = true;
  for (size_t i = 0; i < n_; ++i) {
    Vector2D v1 = (vertices_[i] - vertices_[prev(i)]).normalize();
    Vector2D v2 = (vertices_[next(i)] - vertices_[i]).normalize();
    if ((v2 - v1 * (v1 * v2))[0] * v1[1] > 0)
      convex = false;
  }

  if (!convex)
    return DomainCircular::update();

  return Domain::update();
}

} // namespace Transfinite
