#include <algorithm>

#include "geometry.hh"

namespace Geometry {

BCurve::BCurve() {
}

BCurve::BCurve(const PointVector &cpts) : cp_(cpts) {
}

Point3D
BCurve::eval(double u) const {
  DoubleVector coeff; bernstein(n_, u, coeff);
  Point3D p(0.0, 0.0, 0.0);
  for(size_t k = 0; k <= n_; ++k)
    p += cp_[k] * coeff[k];
  return p;
}

Point3D
BCurve::eval(double u, size_t nr_der, VectorVector &der) const {
  size_t const du = std::min(nr_der, n_);
  der.clear(); der.reserve(nr_der + 1);
  std::vector<DoubleVector> coeff; bernsteinAll(n_, u, coeff);
  std::vector<PointVector> dcp; derivativeControlPoints(du, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector3D(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= n_ - k; ++j)
      der[k] += dcp[k][j] * coeff[n_-k][j];
  }
  for(size_t k = n_ + 1; k <= nr_der; ++k)
    der.push_back(Vector3D(0.0, 0.0, 0.0));
  return der[0];
}

void
BCurve::reverse() {
  std::reverse(cp_.begin(), cp_.end());
}

void
BCurve::normalize() {
}

double
BCurve::arcLength(double from, double to) const {
  if (from >= to)
    return 0.0;

  const static double gauss[] = {-0.861136312, 0.347854845,
                                 -0.339981044, 0.652145155,
                                 0.339981044, 0.652145155,
                                 0.861136312, 0.347854845};

  double sum = 0.0;
  for (size_t i = 0; i < 8; i += 2) {
    VectorVector der;
    double u = ((to - from) * gauss[i] + from + to) * 0.5;
    eval(u, 1, der);
    sum += der[1].norm() * gauss[i+1] * (to - from) * 0.5;
  }

  return sum;
}

void
BCurve::bernstein(size_t n, double u, DoubleVector &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double const u1 = 1.0 - u;
  for(size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for(size_t k = 0; k < j; ++k) {
      double const tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

void
BCurve::bernsteinAll(size_t n, double u, std::vector<DoubleVector> &coeff) {
  coeff.clear(); coeff.resize(n + 1);
  coeff[0].push_back(1.0);
  double const u1 = 1.0 - u;
  for(size_t j = 1; j <= n; ++j) {
    coeff[j].reserve(j + 1);
    double saved = 0.0;
    for(size_t k = 0; k < j; ++k) {
      double const tmp = coeff[j-1][k];
      coeff[j].push_back(saved + tmp * u1);
      saved = tmp * u;
    }
    coeff[j].push_back(saved);
  }
}

void
BCurve::derivativeControlPoints(size_t d, std::vector<PointVector> &dcp) const {
  dcp.clear(); dcp.resize(d + 1);
  dcp[0] = cp_;
  for(size_t k = 1; k <= d; ++k) {
    size_t tmp = n_ - k + 1;
    dcp[k].reserve(tmp);
    for(size_t i = 0; i <= n_ - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp);
  }
}

} // namespace Geometry
