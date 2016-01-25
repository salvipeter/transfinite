#include <array>
#include <cassert>

#include "geometry.hh"

using Geometry::Point2DVector;

struct ControlPointCenters {
  static const Point2DVector &get(size_t n, size_t d) {
    assert(3 <= n && n <= 8);
    assert(3 <= d && d <= 10);
    return data_[n-3][d-3];
  }
private:
  static const Point2DVector data_[6][8];
};
