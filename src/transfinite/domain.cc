#include <algorithm>
#include <cmath>
#include <numeric>

#include "domain.hh"

namespace Transfinite {

Domain::Domain()
  : n_(0) {
}

Domain::~Domain() {
}

void
Domain::setSide(size_t i, const std::shared_ptr<BSCurve> &curve) {
  if (curves_.size() <= i)
    curves_.resize(i + 1);
  curves_[i] = curve;
}

void
Domain::setSides(const CurveVector &curves) {
  curves_ = curves;
}

// Computes everything from vertices
// except for parameters, which are computed only when needed
bool
Domain::update() {
  n_ = vertices_.size();
  computeCenter();
  parameters_.clear();
  du_.resize(n_); dv_.resize(n_);
  for (size_t i = 0; i < n_; ++i) {
    du_[i] = vertices_[i] - vertices_[prev(i)];
    dv_[i] = Vector2D(du_[i][1], -du_[i][0]);
    if ((center_ - vertices_[i]) * dv_[i] < 0)
      dv_[i] = -dv_[i];
  }
  return true;
}

Point2DVector const &
Domain::vertices() const {
  return vertices_;
}

Point2D
Domain::toLocal(size_t i, const Vector2D &v) const {
  return Point2D(v * du_[i], v * dv_[i]) / du_[i].normSqr();
}

const Point2DVector &
Domain::parameters(size_t resolution) const {
  size_t size = 1 + n_ * resolution * (resolution + 1) / 2;
  if (parameters_.size() != size) {
    parameters_.reserve(size);
    parameters_.push_back(center_);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = vertices_[prev(k)] * (1.0 - v) + vertices_[k] * v;
          Point2D p = center_ * (1.0 - u) + ep * u;
          parameters_.push_back(p);
        }
    }
  }
  return parameters_;
}

TriMesh
Domain::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(1 + n_ * resolution * (resolution + 1) / 2);

  size_t inner_start = 0, outer_vert = 1;
  for (size_t layer = 1; layer <= resolution; ++layer) {
    size_t inner_vert = inner_start, outer_start = outer_vert;
    for (size_t side = 0; side < n_; ++side) {
      size_t vert = 0;
      while(true) {
        size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
        mesh.addTriangle(inner_vert, outer_vert, next_vert);
        ++outer_vert;
        if (++vert == layer)
          break;
        size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
        mesh.addTriangle(inner_vert, next_vert, inner_next);
        inner_vert = inner_next;
      }
    }
    inner_start = outer_start;
  }

  return mesh;
}

const Point2D &
Domain::center() const {
  return center_;
}

Point2D
Domain::edgePoint(size_t i, double s) const {
  return vertices_[i] * s + vertices_[prev(i)] * (1.0 - s);
}

double
Domain::edgeLength(size_t i) const {
  return (vertices_[i] - vertices_[prev(i)]).norm();
}

double
Domain::angle(size_t i) const {
  Vector2D v1 = vertices_[prev(i)] - vertices_[i];
  Vector2D v2 = vertices_[next(i)] - vertices_[i];
  return std::acos(std::min(std::max(v1.normalize() * v2.normalize(), -1.0), 1.0));
}

void
Domain::computeCenter() {
  DoubleVector lengths;
  lengths.reserve(n_);
  for (size_t i = 0; i < n_; ++i)
    lengths.push_back((vertices_[i] - vertices_[prev(i)]).norm());
  center_ = Point2D(0, 0);
  for (size_t i = 0; i < n_; ++i)
    center_ += vertices_[i] * (lengths[i] + lengths[next(i)]);
  center_ /= std::accumulate(lengths.begin(), lengths.end(), 0.0) * 2.0;
}

} // namespace Transfinite
