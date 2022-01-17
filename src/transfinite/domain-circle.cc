#include "domain-circle.hh"

namespace Transfinite {

DomainCircle::~DomainCircle() {
}

bool
DomainCircle::update() {
  n_ = curves_.size();
  vertices_.resize(n_);         // other functions depend on this
  return true;
}

size_t
DomainCircle::meshSize(size_t resolution) const {
  return 1 + n_ * resolution * (resolution + 1) / 2;
}

// Similarly to as in Domain, but also for n = 3 & 4
Point2DVector
DomainCircle::parametersImpl(size_t resolution) const {
  size_t size = meshSize(resolution);
  Point2DVector parameters;
  parameters.reserve(size);
  parameters.push_back({0, 0});
  for (size_t j = 1; j <= resolution; ++j) {
    double r = (double)j / (double)resolution;
    for (size_t k = 0; k < n_; ++k)
      for (size_t i = 0; i < j; ++i) {
        double phi = (k + (double)i / j) * 2 * M_PI / n_;
        parameters.push_back({r, phi});
      }
  }
  return parameters;
}

// Same as in Domain, but also for n = 3 & 4
TriMesh
DomainCircle::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(resolution));
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

bool
DomainCircle::onEdge(size_t resolution, size_t index) const {
  size_t size = meshSize(resolution);
  return index >= size - n_ * resolution;
}

void
DomainCircle::computeCenter() {
  center_ = Point2D(0, 0);
}

Point2D
DomainCircle::edgePoint(size_t i, double s) const {
  return { 1, (i + s) * 2 * M_PI / n_ };
}

}
