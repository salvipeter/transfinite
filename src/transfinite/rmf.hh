#pragma once

#include "curves.hh"

namespace Transfinite {

using namespace Geometry;

class RMF {
public:
  void setCurve(const std::shared_ptr<Curve> &c);
  void setStart(const Vector3D &start);
  void setEnd(const Vector3D &end);
  void update();
  Vector3D eval(double u) const;

private:
  struct Frame {
    Frame(double u, double s, const Point3D &p, const Vector3D &d, const Vector3D &n)
      : u(u), s(s), p(p), d(d), n(n) {
    }
    double u, s;
    Point3D p;
    Vector3D d, n;
  };
  using Matrix3x3 = std::array<std::array<double, 3>, 3>;
  static const Matrix3x3 &rotationMatrix(const Vector3D &u, double theta);
  static void rotateFrame(Frame &f, double angle);
  Frame nextFrame(const Frame &prev, double u) const;

  const static size_t resolution_;
  std::shared_ptr<Curve> curve_;
  Vector3D start_, end_;
  std::vector<Frame> frames_;
  double angleCorrection_;
};

} // namespace Transfinite
