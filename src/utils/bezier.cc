#include "bezier.hh"

#include <algorithm>
#include <cmath>

#include "nelder-mead.hh"

namespace Geometry {

BCurve::BCurve() {
}

BCurve::BCurve(const PointVector &cpts) : n_(cpts.size() - 1), cp_(cpts)
{
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

const PointVector &
BCurve::controlPoints() const {
  return cp_;
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
BCurve::fitClassA(size_t degree,
                  const Point3D &pa, const Vector3D &va,
                  const Point3D &pb, const Vector3D &vb) {
  // Initialize values
  n_ = degree;
  cp_.resize(degree + 1);
  Vector3D v0 = va; v0.normalize();
  Vector3D v1 = vb; v1.normalize();

  // Generating control points
  auto generateClassA = [this] (const Point3D &p, const Vector3D &v, const Matrix3x3 &M) {
    cp_[0] = p;
    cp_[1] = p + v;
    Vector3D next = v;
    for (size_t i = 2; i <= n_; ++i) {
      next = M * next;
      cp_[i] = cp_[i-1] + next;
    }
  };
  auto setupMatrix = [v0, v1, degree] (double phi, double s, Matrix3x3 &M) {
    // TODO: corner case - when v0 and v1 are parallel directions
    Vector3D axis = (Matrix3x3::rotation((v1 - v0).normalize(), phi) * (v0 ^ v1)).normalize();
    Vector3D w0 = (v0 - axis * (v0 * axis)).normalize();
    Vector3D w1 = (v1 - axis * (v1 * axis)).normalize();
    double theta = std::acos(std::min(std::max(w0 * w1, -1.0), 1.0)) / (degree - 1);
    M = Matrix3x3::rotation(axis, theta) * s;
  };
  auto evalMatrix = [this, generateClassA, pa, pb, v0] (const Matrix3x3 &M) -> double {
    generateClassA(pa, v0, M);
    Vector3D u1 = (cp_[n_] - cp_[0]).normalize();
    Vector3D u2 = (pb - pa).normalize();
    return u1 * u2 - 1.0;
  };

  // Optimization
  Matrix3x3 M;
  auto targetFunction = [setupMatrix, evalMatrix, &M] (const std::vector<double> &phi_s) -> double {
    setupMatrix(phi_s[0], phi_s[1], M);
    return std::abs(evalMatrix(M));
  };
  std::vector<double> start = { 0, 1 };
  NelderMead::optimize(targetFunction, start, 100, 1.0e-15, 1.0);
  double err = targetFunction(start);        // update to the best value found so far
  if (err > 0.001 && degree < 50) {
    fitClassA(degree * 2, pa, va, pb, vb);
    return;
  }

  // Postprocessing
  v0 *= (pb - pa).norm() / (cp_[n_] - cp_[0]).norm();
  generateClassA(pa, v0, M);
  cp_[degree] = pb;
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
