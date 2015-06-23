#pragma once

#include <memory>

class Vector2D {
public:
  // Constructors & destructor
  Vector2D();
  Vector2D(double x, double y);
  Vector2D(const Vector2D &v);
  Vector2D(Vector2D &&v);
  ~Vector2D();

  // Assignments
  Vector2D &operator=(const Vector2D &v);
  Vector2D &operator=(Vector2D &&v);
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
  class Vector2DImpl;
  Vector2D(std::unique_ptr<Vector2DImpl> &&impl);
  std::unique_ptr<Vector2DImpl> impl_;
};

using Point2D = Vector2D;

class Vector3D {
public:
  // Constructors & destructor
  Vector3D();
  Vector3D(double x, double y, double z);
  Vector3D(const Vector3D &v);
  Vector3D(Vector3D &&v);
  ~Vector3D();

  // Assignments
  Vector3D &operator=(const Vector3D &v);
  Vector3D &operator=(Vector3D &&v);
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
  class Vector3DImpl;
  Vector3D(std::unique_ptr<Vector3DImpl> &&impl);
  std::unique_ptr<Vector3DImpl> impl_;
};

using Point3D = Vector3D;

// class BSCurve {
// public:
//   BSplineCurve<Point3D> BSCurve;
// };
