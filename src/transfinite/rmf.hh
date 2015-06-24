#pragma once

#include <array>

#include "geometry.hh"

class RMF
{
public:
  RMF();
  RMF(const BSCurve &c, const Vector3D &start, const Vector3D &end);
  Vector3D eval(double u) const;
private:
  struct Frame {
    Frame(double u, double s, const Point3D &p, const Vector3D &d, const Vector3D &n)
      : u(u), s(s), p(p), d(d), n(n) {}
    double u, s;
    Point3D p;
    Vector3D d, n;
  };
  using Matrix3x3 = std::array<std::array<double, 3>, 3>;
  static const Matrix3x3 &rotationMatrix(const Vector3D &u, double theta);
  static void rotateFrame(Frame &f, double angle);
  static Frame nextFrame(const BSCurve &c, const Frame &prev, double u);

  const static size_t resolution_;
  BSCurve curve_;
  std::vector<Frame> frames_;
  double angleCorrection_;
};
