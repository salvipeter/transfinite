#include <algorithm>
#include <cmath>

#include "rmf.hh"
#include "utilities.hh"

namespace Transfinite {

const size_t RMF::resolution_ = 100;

void
RMF::setCurve(const std::shared_ptr<Curve> &c) {
  curve_ = c;
}

void
RMF::setStart(const Vector3D &start) {
  start_ = start;
}

void
RMF::setEnd(const Vector3D &end) {
  end_ = end;
}

void
RMF::update() {
  // As described in `Computation of Rotation Minimizing Frames', Wang et al., 2008.
  // A limitation of this method is that it is determined by the starting frame
  // and the curve tangents, so an end frame cannot be supplied.
  frames_.clear();
  frames_.reserve(resolution_ + 1);
  VectorVector der; curve_->eval(0.0, 1, der);
  Frame f(0.0, 0.0, der[0], der[1].normalize(), start_);
  for (size_t i = 1; i <= resolution_; ++i) {
    double u = (double)i / (double)resolution_;
    frames_.push_back(f);
    f = nextFrame(f, u);
  }
  frames_.push_back(f);

  // As a workaround, we can add a rotation gradually,
  // by minimizing the total squared angular speed, see section 6.3 in the paper.
  const Vector3D &rmfEnd = f.n;
  angleCorrection_ = std::acos(inrange(-1, end_ * rmfEnd, 1));
  if (((rmfEnd - end_) ^ end_) * f.d < 0.0)
    angleCorrection_ *= -1.0;
  angleCorrection_ /= f.s;
}

Vector3D
RMF::eval(double u) const {
  auto i = std::upper_bound(frames_.begin(), frames_.end(), u,
                            [](double x, const Frame &f) { return x < f.u; });
  Frame f = nextFrame(*(--i), u);
  rotateFrame(f, f.s * angleCorrection_);
  return f.n;
}

const RMF::Matrix3x3 &
RMF::rotationMatrix(const Vector3D &u, double theta) {
  static Matrix3x3 m;

  double x = u[0], y = u[1], z = u[2];
  double c = std::cos(theta), c1 = 1.0 - c, s = std::sin(theta);
  m[0][0] = c + x * x * c1;
  m[1][0] = y * x * c1 + z * s;
  m[2][0] = z * x * c1 - y * s;
  m[0][1] = x * y * c1 - z * s;
  m[1][1] = c + y * y * c1;
  m[2][1] = z * y * c1 + x * s;
  m[0][2] = x * z * c1 + y * s;
  m[1][2] = y * z * c1 - x * s;
  m[2][2] = c + z * z * c1;

  return m;
}

void
RMF::rotateFrame(Frame &f, double angle) {
  const Matrix3x3 &r = rotationMatrix(f.d, angle);
  Vector3D n(0.0, 0.0, 0.0);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      n[i] += r[i][j] * f.n[j];
  f.n = n;
}

RMF::Frame
RMF::nextFrame(const Frame &prev, double u) const {
  VectorVector der; curve_->eval(u, 1, der);
  der[1].normalize();
  Vector3D v1 = der[0] - prev.p;
  double c1 = v1.normSqr();
  if (c1 < epsilon)
    return prev;
  Vector3D v2 = v1 * 2 / c1;
  Vector3D nL = prev.n - v2 * (v1 * prev.n);
  Vector3D dL = prev.d - v2 * (v1 * prev.d);
  v2 = der[1] - dL;
  double c2 = v2.normSqr();
  Vector3D nNext;
  if (c2 < epsilon)
    nNext = nL;
  else
    nNext = nL - v2 * 2 / c2 * (v2 * nL);
  return Frame(u, prev.s + curve_->arcLength(prev.u, u), der[0], der[1], nNext);
}

} // namespace Transfinite
