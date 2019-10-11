#include <cmath>

#include "utilities.hh"

namespace Transfinite {

Point2D
affineCombine(const Point2D &p, double x, const Point2D &q) {
  return p * (1 - x) + q * x;
}

Point3D
affineCombine(const Point3D &p, double x, const Point3D &q) {
  return p * (1 - x) + q * x;
}

double
inrange(double min, double x, double max) {
  if (x < min)
    return min;
  if (x > max)
    return max;
  return x;
}

size_t
binomial(size_t n, size_t k) {
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

double
hermite(int i, double t) {
  switch(i) {
  case 0: return std::pow(1 - t, 3) + 3.0 * std::pow(1 - t, 2) * t;
  case 1: return std::pow(1 - t, 2) * t;
  case 2: return (1 - t) * std::pow(t, 2);
  case 3: return 3.0 * (1 - t) * std::pow(t, 2) + std::pow(t, 3);
  default: ;
  }
  return -1.0;                  // should not come here
}

void
bernstein(size_t n, double u, DoubleVector &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double  tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

double
bernstein(size_t i, size_t n, double u) {
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  const double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t j = n; j >= k; --j)
      tmp[j] = tmp[j-1] * u + tmp[j] * u1;
  return tmp[n];
}

void
bezierElevate(PointVector &cpts) {
  size_t n = cpts.size();
  Point3D tmp = cpts[0];
  for (size_t i = 1; i < n; ++i) {
    tmp = tmp * i / n + cpts[i] * (n - i) / n;
    std::swap(cpts[i], tmp);
  }
  cpts.push_back(tmp);
}

} // namespace Transfinite
