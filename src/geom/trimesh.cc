#include "geometry.hh"

#include <iostream>
#include <fstream>
#include <string>

namespace Geometry {

void
TriMesh::resizePoints(size_t n) {
  points_.resize(n);
}

void
TriMesh::setPoint(size_t i, const Point3D &p) {
  points_[i] = p;
}

void
TriMesh::setPoints(const PointVector &pv) {
  points_ = pv;
}

void
TriMesh::addTriangle(size_t a, size_t b, size_t c) {
  triangles_.push_back({a, b, c});
}

PointVector
TriMesh::points() const {
  return points_;
}

std::list<TriMesh::Triangle>
TriMesh::triangles() const {
  return triangles_;
}

void
TriMesh::writeOBJ(std::string filename) const {
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &p : points_)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (const auto &t : triangles_)
    f << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << std::endl;
  f.close();
}

} // namespace Geometry
