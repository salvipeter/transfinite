#include <cmath>

#include "domain.hh"

Domain::Domain() : n_(0) {
}

Domain::~Domain() {
}

// Computes everything from vertices
// except for parameters, which are computed only when needed
void
Domain::invalidate() {
  n_ = vertices_.size();
  computeCenter();
  parameters_.clear();
  du_.resize(n_); dv_.resize(n_);
  for(size_t i = 0; i < n_; ++i) {
    du_[i] = vertices_[i] - vertices_[prev(i)];
    dv_[i] = Vector2D(du_[i][1], -du_[i][0]);
    if((center_ - vertices_[i]) * dv_[i] < 0)
      dv_[i] = -dv_[i];
  }
  local_vertices_.resize(n_);
  for(size_t i = 0; i < n_; ++i) {
    local_vertices_[i].resize(n_);
    for(size_t j = 0; j < n_; ++j)
      local_vertices_[i][j] = toLocal(i, vertices_[j]);
  }
}

Point2DVector const &
Domain::verticesGlobal() const {
  return vertices_;
}

Point2DVector const &
Domain::verticesLocal(size_t i) const {
  return local_vertices_[i];
}

Point2D
Domain::toLocal(size_t i, const Point2D &p) const {
  Vector2D const v = p - vertices_[prev(i)];
  if(std::abs(dv_[i][1]) < epsilon) {
    double const a = (v[1] - v[0] * dv_[i][1] / dv_[i][0]) /
      (du_[i][1] - du_[i][0] * dv_[i][1] / dv_[i][0]);
    return Point2D(v[0] - a * du_[i][0] - dv_[i][0], a);
  }
  double const a = (v[0] - v[1] * dv_[i][0] / dv_[i][1]) /
    (du_[i][0] - du_[i][1] * dv_[i][0] / dv_[i][1]);
  return Point2D(a, v[1] - a * du_[i][1] - dv_[i][1]);
}

Point2D
Domain::toGlobal(size_t i, const Point2D &p) const {
  return vertices_[prev(i)] + du_[i] * p[0] + dv_[i] * p[1];
}

const Point2DVector &
Domain::globalParameters(size_t resolution) const {
  size_t size = 1 + n_ * resolution * (resolution + 1) / 2;
  if(parameters_.size() != size) {
    parameters_.reserve(size);
    parameters_.push_back(center_);
    for(size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for(size_t k = 0; k < n_; ++k)
        for(size_t i = 0; i < j; ++i) {
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
  for(size_t layer = 1; layer <= resolution; ++layer) {
    size_t inner_vert = inner_start, outer_start = outer_vert;
    for(size_t side = 0; side < n_; ++side) {
      size_t vert = 0;
      while(true) {
        size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
        mesh.addTriangle(inner_vert, outer_vert, next_vert);
        ++outer_vert;
        if(++vert == layer)
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
