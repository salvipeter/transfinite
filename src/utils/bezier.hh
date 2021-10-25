#pragma once

#include <geometry.hh>

namespace Geometry {

class BCurve {
public:
  // Constructors
  BCurve();
  BCurve(const PointVector &cpts);

  // Evaluation
  Point3D eval(double u) const;
  Point3D eval(double u, size_t nr_der, VectorVector &der) const;

  // Coordinates
  const PointVector &controlPoints() const;

  // Parameterization
  void reverse();
  void normalize();

  // Other
  void fitClassA(size_t degree,
                 const Point3D &pa, const Vector3D &va,
                 const Point3D &pb, const Vector3D &vb);
  double arcLength(double from, double to) const;

private:
  static void bernstein(size_t n, double u, DoubleVector &coeff);
  static void bernsteinAll(size_t n, double u, std::vector<DoubleVector> &coeff);
  void derivativeControlPoints(size_t d, std::vector<PointVector> &dcp) const;

  size_t n_;
  PointVector cp_;
};

}
