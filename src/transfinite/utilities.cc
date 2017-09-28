#include <cmath>

#include "utilities.hh"

namespace Transfinite {

double
blendHermite(double x) {
  double x2 = x * x;
  return 2.0 * x * x2 - 3.0 * x2 + 1.0;
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

} // namespace Transfinite
