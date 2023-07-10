#include <algorithm>
#include <cmath>
#include <numeric>

#include "domain.hh"
#include "utilities.hh"

namespace Transfinite {

Domain::Domain()
  : n_(0) {
}

Domain::~Domain() {
}

void
Domain::setSide(size_t i, const std::shared_ptr<Curve> &curve) {
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

Point2D
Domain::fromLocal(size_t i, const Vector2D &v) const {
  return du_[i] * v[0] + dv_[i] * v[1];
}

bool
Domain::intersectEdgeWithRay(size_t i, const Point2D &p, const Vector2D &v, Point2D &result) const {
  Point2D q1 = vertices_[prev(i)], q2 = vertices_[i];
  double denom = v[0] * (q2[1] - q1[1]) - v[1] * (q2[0] - q1[0]);

  if (fabs(denom) < epsilon)
    return false;

  double t = (v[0] * (p[1] - q1[1]) - v[1] * (p[0] - q1[0])) / denom;

  if (t < -epsilon || t > 1 + epsilon)
    return false;

  t = inrange(0.0, t, 1.0);
  result = q1 * (1 - t) + q2 * t;

  return true;
}

size_t
Domain::size() const {
  return n_;
}

size_t
Domain::meshSize(size_t resolution) const {
  if (n_ == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n_ == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n_ * resolution * (resolution + 1) / 2;
}

const Point2DVector &
Domain::parameters(size_t resolution) const {
  size_t size = meshSize(resolution);
  if (parameters_.size() == size)
    return parameters_;
  parameters_ = parametersImpl(resolution);
  return parameters_;
}

Point2DVector
Domain::parametersImpl(size_t resolution) const {
  Point2DVector parameters;
  parameters.reserve(meshSize(resolution));

  if (n_ == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = vertices_[0] * u + vertices_[2] * (1 - u);
      auto q = vertices_[1] * u + vertices_[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        parameters.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n_ == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = vertices_[0] * (1 - u) + vertices_[1] * u;
      auto q = vertices_[3] * (1 - u) + vertices_[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        parameters.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n_ > 4
    parameters.push_back(center_);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = vertices_[prev(k)] * (1.0 - v) + vertices_[k] * v;
          Point2D p = center_ * (1.0 - u) + ep * u;
          parameters.push_back(p);
        }
    }
  }
  return parameters;
}

bool
Domain::onEdge(size_t resolution, size_t index) const {
  if (n_ == 3) {
    if (index >= meshSize(resolution) - resolution - 1)
      return true;
    auto issquare = [](size_t n) {
                      size_t root = std::round(std::sqrt(n));
                      return root * root == n;
                    };
    size_t n = index * 8 + 1;
    return issquare(n) || issquare(n + 8);
  }
  if (n_ == 4) {
    return index <= resolution || index >= (resolution + 1) * resolution ||
      index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
  }
  return index >= meshSize(resolution) - n_ * resolution;
}

TriMesh
Domain::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(resolution));

  if (n_ == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n_ == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n_ > 4
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
  return std::acos(inrange(-1, v1.normalize() * v2.normalize(), 1));
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
