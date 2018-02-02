#include "geometry.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace Geometry {

void
TriMesh::clear() {
  points_.clear();
  triangles_.clear();
}

void
TriMesh::resizePoints(size_t n) {
  points_.resize(n);
}

void
TriMesh::setPoints(const PointVector &pv) {
  points_ = pv;
}

void
TriMesh::addTriangle(size_t a, size_t b, size_t c) {
  triangles_.push_back({a, b, c});
}

void
TriMesh::setTriangles(const std::list<Triangle> &tl) {
  triangles_ = tl;
}

Point3D &
TriMesh::operator[](size_t i) {
  return points_[i];
}
 
const Point3D &
TriMesh::operator[](size_t i) const {
  return points_[i];
}

PointVector
TriMesh::points() const {
  return points_;
}

std::list<TriMesh::Triangle>
TriMesh::triangles() const {
  return triangles_;
}

Point3D
TriMesh::projectToTriangle(const Point3D &p, const Triangle &tri) const {
  const Point3D &q1 = points_[tri[0]], &q2 = points_[tri[1]], &q3 = points_[tri[2]];
  // As in Schneider, Eberly: Geometric Tools for Computer Graphics, Morgan Kaufmann, 2003.
  // Section 10.3.2, pp. 376-382 (with my corrections)
  const Point3D &P = p, &B = q1;
  Vector3D E0 = q2 - q1, E1 = q3 - q1, D = B - P;
  double a = E0 * E0, b = E0 * E1, c = E1 * E1, d = E0 * D, e = E1 * D;
  double det = a * c - b * b, s = b * e - c * d, t = b * d - a * e;
  if (s + t <= det) {
    if (s < 0) {
      if (t < 0) {
        // Region 4
        if (e < 0) {
          s = 0.0;
          t = (-e >= c ? 1.0 : -e / c);
        } else if (d < 0) {
          t = 0.0;
          s = (-d >= a ? 1.0 : -d / a);
        } else {
          s = 0.0;
          t = 0.0;
        }
      } else {
        // Region 3
        s = 0.0;
        t = (e >= 0.0 ? 0.0 : (-e >= c ? 1.0 : -e / c));
      }
    } else if (t < 0) {
      // Region 5
      t = 0.0;
      s = (d >= 0.0 ? 0.0 : (-d >= a ? 1.0 : -d / a));
    } else {
      // Region 0
      double invDet = 1.0 / fabs(det);
      s *= invDet;
      t *= invDet;
    }
  } else {
    if (s < 0) {
      // Region 2
      double tmp0 = b + d, tmp1 = c + e;
      if (tmp1 > tmp0) {
        double numer = tmp1 - tmp0;
        double denom = a - 2 * b + c;
        s = (numer >= denom ? 1.0 : numer / denom);
        t = 1.0 - s;
      } else {
        s = 0.0;
        t = (tmp1 <= 0.0 ? 1.0 : (e >= 0.0 ? 0.0 : -e / c));
      }
    } else if (t < 0) {
      // Region 6
      double numer = c + e - b - d, denom = c - 2 * b + a;
      if (numer < denom) {
        s = numer / denom;
        t = 1.0 - s;
      } else {
        t = 0.0;
        s = (d >= 0.0 ? 0.0 : (-d >= a ? 1.0 : -d / a));
      }
    } else {
      // Region 1
      double numer = c + e - b - d;
      if (numer <= 0) {
        s = 0;
      } else {
        double denom = a - 2 * b + c;
        s = (numer >= denom ? 1.0 : numer / denom);
      }
      t  = 1.0 - s;
    }
  }
  return B + E0 * s + E1 * t;
}

TriMesh::Triangle
TriMesh::closestTriangle(const Point3D &p) const {
  // Trivial (slow) implementation
  std::list<Triangle>::const_iterator i = triangles_.begin(), result = i;
  double min = (projectToTriangle(p, *i) - p).norm();
  while (++i != triangles_.end()) {
    double d = (projectToTriangle(p, *i) - p).norm();
    if (d < min) {
      min = d;
      result = i;
    }
  }
  return *result;
}

bool
TriMesh::writeOBJ(std::string filename) const {
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return false;
  }
  for (const auto &p : points_)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (const auto &t : triangles_)
    f << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << std::endl;
  f.close();
  return true;
}

} // namespace Geometry
