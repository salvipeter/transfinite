#pragma once

#include "geometry.hh"

#include "Vector.hh"

class Vector2D::Vector2DImpl {
public:
  Vector2DImpl() {}
  Vector2DImpl(double x, double y) : v_(Vector<2, double>(x, y)) {}
  Vector2DImpl &operator+=(Vector2DImpl v) { v_ += v.v_; return *this; }
  Vector2DImpl &operator-=(Vector2DImpl v) { v_ -= v.v_; return *this; }
  Vector2DImpl &operator*=(double x) { v_ *= x; return *this; }
  Vector2DImpl &operator/=(double x) { v_ /= x; return *this; }
  double &operator[](size_t i) { return v_[i]; }
  const double &operator[](size_t i) const { return v_[i]; }
  Vector2DImpl operator-() const { return Vector2DImpl(-v_); }
  Vector2DImpl operator+(const Vector2DImpl &v) const { return Vector2DImpl(v_ + v.v_); }
  Vector2DImpl operator-(const Vector2DImpl &v) const { return Vector2DImpl(v_ - v.v_); }
  double operator*(const Vector2DImpl &v) const { return v_ * v.v_; }
  Vector2DImpl operator*(double x) const { return Vector2DImpl(v_ * x); }
  Vector2DImpl operator/(double x) const { return Vector2DImpl(v_ / x); }
  double norm() const { return v_.Norm(); }
  double normSqr() const { return v_.NormSqr(); }
  Vector2DImpl &normalize() { v_.Normalize(); return *this; }
private:
  Vector2DImpl(Vector<2, double> v) : v_(v) {}
  Vector<2, double> v_;
};

class Vector3D::Vector3DImpl {
public:
  Vector3DImpl() {}
  Vector3DImpl(double x, double y, double z) : v_(Vector<3, double>(x, y, z)) {}
  Vector3DImpl &operator+=(Vector3DImpl v) { v_ += v.v_; return *this; }
  Vector3DImpl &operator-=(Vector3DImpl v) { v_ -= v.v_; return *this; }
  Vector3DImpl &operator*=(double x) { v_ *= x; return *this; }
  Vector3DImpl &operator/=(double x) { v_ /= x; return *this; }
  double &operator[](size_t i) { return v_[i]; }
  const double &operator[](size_t i) const { return v_[i]; }
  Vector3DImpl operator-() const { return Vector3DImpl(-v_); }
  Vector3DImpl operator+(const Vector3DImpl &v) const { return Vector3DImpl(v_ + v.v_); }
  Vector3DImpl operator-(const Vector3DImpl &v) const { return Vector3DImpl(v_ - v.v_); }
  Vector3DImpl operator^(const Vector3DImpl &v) const { return Vector3DImpl(v_ ^ v.v_); }
  double operator*(const Vector3DImpl &v) const { return v_ * v.v_; }
  Vector3DImpl operator*(double x) const { return Vector3DImpl(v_ * x); }
  Vector3DImpl operator/(double x) const { return Vector3DImpl(v_ / x); }
  double norm() const { return v_.Norm(); }
  double normSqr() const { return v_.NormSqr(); }
  Vector3DImpl &normalize() { v_.Normalize(); return *this; }
private:
  Vector3DImpl(Vector<3, double> v) : v_(v) {}
  Vector<3, double> v_;
};
