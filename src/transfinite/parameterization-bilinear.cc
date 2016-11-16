#include <cmath>

#include "domain.hh"
#include "parameterization-bilinear.hh"

namespace Transfinite {

ParameterizationBilinear::~ParameterizationBilinear() {
}

Point2D
ParameterizationBilinear::mapToRibbon(size_t i, const Point2D &uv) const {
  Point2DVector const &vertices = domain_->vertices();

  const Point2D &v0 = vertices[prev(i, 2)];
  const Point2D &v1 = vertices[prev(i)];
  const Point2D &v2 = vertices[i];
  const Point2D &v3 = vertices[next(i)];

  // bilinear: p1 - (0,0) - (1,0) - p2
  Point2D p1 = domain_->toLocal(i, v0 - v1);
  Point2D p2 = domain_->toLocal(i, v3 - v2);
  Point2D p = domain_->toLocal(i, uv - v1);

  double a = p1[1] - p2[1];
  double b = p[0] * p2[1] - (p[0] + 1.0) * p1[1] + (p1[0] - p2[0]) * p[1];
  double c = p[0] * p1[1] - p[1] * p1[0];

  Point2D sd;

  // solve: a s^2 + b s + c = 0
  if (std::abs(a) < epsilon) {
    if (std::abs(b) < epsilon) {
      // should not come here
      sd[0] = 0.0;
    } else
      sd[0] = -c / b;
  } else {
    double D = b * b - 4.0 * a * c;
    if (D < epsilon)
      D = 0.0;
    D = std::sqrt(D) / (2.0 * a);
    double s = -b / (2.0 * a);

    double s1 = s + D;
    double s2 = s - D;
    if (s1 < s2)
      std::swap(s1, s2);

    double dev1 = std::max(0.0, s1 - 1.0);
    double dev2 = std::max(0.0, -s2);
    if (dev1 < dev2)
      sd[0] = s1;
    else
      sd[0] = s2;
  }

  sd[1] = p[1] / (p1[1] * (1.0 - sd[0]) + p2[1] * sd[0]);
  return sd;
}

} // namespace Transfinite
