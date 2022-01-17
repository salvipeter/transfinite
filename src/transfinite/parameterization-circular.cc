#include "parameterization-circular.hh"

#include <cmath>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Transfinite {

ParameterizationCircular::~ParameterizationCircular() {
}

static std::pair<Point2D, double> hcircle(size_t n, double h) {
  double phi = (2 * h + 1) * M_PI / n;
  Point2D p(cos(phi), sin(phi));
  Vector2D t(sin(h * M_PI - phi), cos(h * M_PI - phi));
  if (std::abs(t[0]) < epsilon)
    return { p, 0 };            // line segment
  Point2D O(p[0] + p[1] * t[1] / t[0], 0);
  return { O, (O - p).norm() };
}

static double bisect(const std::function<double(double)> &f, double xl, double xh) {
  double xm = xl;
  for (size_t i = 0; i < 100; ++i) {
    double xm_old = xm;
    xm = (xl + xh) / 2;
    if (xm != 0 && std::abs((xm - xm_old) / xm) < epsilon)
      break;
    double test = f(xl) * f(xm);
    if (test < 0)
      xh = xm;
    else if (test > 0)
      xl = xm;
    else
      break;
  }
  return xm;
}

Point2D
ParameterizationCircular::mapToRibbon(size_t i, const Point2D &uv) const {
  // uv is actually (r, phi), with the first segment between phi = 0 and phi = 2pi/n
  // Here we first compute a local phi by rotating the point
  double r = uv[0];
  double phi = std::fmod(uv[1] + (3 - (i * 2 + 1.0) / n_) * M_PI, 2 * M_PI) - M_PI;
  if (r > 1 - epsilon) {
    if (std::abs(phi) <= M_PI / n_)
      return { 0, 0 };
    if (std::abs(phi) >= 3 * M_PI / n_)
      return { 1, 1 };
    double h = (std::abs(phi) * n_ / M_PI - 1) / 2;
    return { h, h };
  }
  Point2D p(r * cos(phi), r * sin(phi));
  auto deviation = [&](double h) {
    auto [O, R] = hcircle(n_, h);
    if (R == 0) // line segment
      return O[0] - p[0];
    if (h > 0.8)                // hack
      return (O - p).norm() - R;
    double D = std::sqrt(R * R - (p[1] - O[1]) * (p[1] - O[1]));
    double x = O[0] + D;
    if (x > 1)
      x = O[0] - D;
    return x - p[0];
  };
  double h;
  if (p[0] < 0)
    h = bisect(deviation, 0.2, 1);
  else
    h = bisect(deviation, 0, 0.8);
  return { h, h };
}

}
