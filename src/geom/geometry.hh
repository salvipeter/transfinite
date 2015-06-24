#pragma once

#include <array>
#include <memory>
#include <vector>

using DoubleVector = std::vector<double>;

class Vector2D {
public:
  // Constructors
  Vector2D();
  Vector2D(double x, double y);

  // Assignments
  Vector2D &operator+=(const Vector2D &v);
  Vector2D &operator-=(const Vector2D &v);
  Vector2D &operator*=(double x);
  Vector2D &operator/=(double x);

  // Coordinates
  double &operator[](size_t i);
  const double &operator[](size_t i) const;

  // Arithmetic
  Vector2D operator-() const;
  Vector2D operator+(const Vector2D &v) const;
  Vector2D operator-(const Vector2D &v) const;
  double operator*(const Vector2D &v) const;
  Vector2D operator*(double x) const;
  Vector2D operator/(double x) const;

  // Other
  double norm() const;
  double normSqr() const;
  Vector2D &normalize();

private:
  std::array<double, 2> v_;
};

using Point2D = Vector2D;

class Vector3D {
public:
  // Constructors
  Vector3D();
  Vector3D(double x, double y, double z);

  // Assignments
  Vector3D &operator+=(const Vector3D &v);
  Vector3D &operator-=(const Vector3D &v);
  Vector3D &operator*=(double x);
  Vector3D &operator/=(double x);

  // Coordinates
  double &operator[](size_t i);
  const double &operator[](size_t i) const;

  // Arithmetic
  Vector3D operator-() const;
  Vector3D operator+(const Vector3D &v) const;
  Vector3D operator-(const Vector3D &v) const;
  Vector3D operator^(const Vector3D &v) const;
  double operator*(const Vector3D &v) const;
  Vector3D operator*(double x) const;
  Vector3D operator/(double x) const;

  // Other
  double norm() const;
  double normSqr() const;
  Vector3D &normalize();

private:
  std::array<double, 3> v_;
};

using VectorVector = std::vector<Vector3D>;
using Point3D = Vector3D;
using PointVector = std::vector<Point3D>;

class BSCurve {
public:
  // Constructors
  BSCurve();
  BSCurve(size_t degree, DoubleVector knots, PointVector cpts);

  // Evaluation
  Point3D eval(double u) const;
  Point3D eval(double u, size_t nr_der, VectorVector &der) const;

  // Parameterization
  void reverse();
  void normalize();

private:
  size_t findSpan(double u) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  void basisFunctionsAll(size_t i, double u, std::vector<DoubleVector> &coeff) const;
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, std::vector<PointVector> &dcp) const;

  size_t p_;
  size_t n_;
  DoubleVector knots_;
  PointVector cp_;
};
