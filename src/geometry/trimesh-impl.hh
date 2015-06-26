#pragma once

#include <algorithm>
#include <fstream>

#include "geometry.hh"

#include "Triangulation.hh"

class TriMesh::TriMeshImpl {
public:
  void resizePoints(size_t n) { m_.Points().resize(n); }
  void setPoint(size_t i, const Point3D &p) { m_.Points()[i] = Point<3, double>(p[0], p[1], p[2]); }
  void setPoints(const PointVector &pv) {
    m_.Points().clear(); m_.Points().reserve(pv.size());
    std::transform(pv.begin(), pv.end(), std::back_inserter(m_.Points()),
                   [](const Point3D &p) { return Point<3, double>(p[0], p[1], p[2]); });
  }
  void addTriangle(size_t a, size_t b, size_t c) { m_.AddTriangle(a, b, c, false); }
  void writeOBJ(std::string filename) const {
    std::ofstream f(filename);
    if (!f.is_open()) {
      std::cerr << "Unable to open file: " << filename << std::endl;
      return;
    }
    int vertex_offs = 0;
    m_.TriangulationToOBJ(f, 0, vertex_offs);
    f.close();
  }
private:
  Triangulation<Point<3, double>> m_;
};
