#include <limits>

#include "vector-impl.hh"
#include "bspline-impl.hh"
#include "trimesh-impl.hh"

//const double EPS_REAL = std::numeric_limits<double>::epsilon();

Vector2D::Vector2D()
  : impl_(std::make_unique<Vector2DImpl>()) {
}

Vector2D::Vector2D(double x, double y)
  : impl_(std::make_unique<Vector2DImpl>(x, y)) {
}

Vector2D::Vector2D(std::unique_ptr<Vector2DImpl> &&impl)
  : impl_(std::move(impl)) {
}

Vector2D::Vector2D(const Vector2D &v)
  : impl_(std::make_unique<Vector2DImpl>(*v.impl_)) {
}

Vector2D::Vector2D(Vector2D &&v) = default;

Vector2D::~Vector2D() = default;

Vector2D &
Vector2D::operator=(Vector2D &&v) = default;

Vector2D &
Vector2D::operator=(const Vector2D &v) {
  *impl_ = *v.impl_;
  return *this;
}

Vector2D &
Vector2D::operator+=(const Vector2D &v) {
  *impl_ += *v.impl_;
  return *this;
}

Vector2D &
Vector2D::operator-=(const Vector2D &v) {
  *impl_ -= *v.impl_;
  return *this;
}

Vector2D &
Vector2D::operator*=(double x) {
  *impl_ *= x;
  return *this;
}

Vector2D &
Vector2D::operator/=(double x) {
  *impl_ /= x;
  return *this;
}

double &
Vector2D::operator[](size_t i) {
  return (*impl_)[i];
}

const double &
Vector2D::operator[](size_t i) const {
  return (*impl_)[i];
}

Vector2D
Vector2D::operator-() const {
  return Vector2D(std::make_unique<Vector2DImpl>(-(*impl_)));
}

Vector2D
Vector2D::operator+(const Vector2D &v) const {
  return Vector2D(std::make_unique<Vector2DImpl>(*impl_ + *v.impl_));
}

Vector2D
Vector2D::operator-(const Vector2D &v) const {
  return Vector2D(std::make_unique<Vector2DImpl>(*impl_ - *v.impl_));
}

double
Vector2D::operator*(const Vector2D &v) const {
  return *impl_ * *v.impl_;
}

Vector2D
Vector2D::operator*(double x) const {
  return Vector2D(std::make_unique<Vector2DImpl>(*impl_ * x));
}

Vector2D
Vector2D::operator/(double x) const {
  return Vector2D(std::make_unique<Vector2DImpl>(*impl_ / x));
}

double
Vector2D::norm() const {
  return impl_->norm();
}

double
Vector2D::normSqr() const {
  return impl_->normSqr();
}

Vector2D &
Vector2D::normalize() {
  impl_->normalize();
  return *this;
}

Vector3D::Vector3D()
  : impl_(std::make_unique<Vector3DImpl>()) {
}

Vector3D::Vector3D(double x, double y, double z)
  : impl_(std::make_unique<Vector3DImpl>(x, y, z)) {
}

Vector3D::Vector3D(std::unique_ptr<Vector3DImpl> &&impl)
  : impl_(std::move(impl)) {
}

Vector3D::Vector3D(const Vector3D &v)
  : impl_(std::make_unique<Vector3DImpl>(*v.impl_)) {
}

Vector3D::Vector3D(Vector3D &&v) = default;

Vector3D::~Vector3D() = default;

Vector3D &
Vector3D::operator=(Vector3D &&v) = default;

Vector3D &
Vector3D::operator=(const Vector3D &v) {
  *impl_ = *v.impl_;
  return *this;
}

Vector3D &
Vector3D::operator+=(const Vector3D &v) {
  *impl_ += *v.impl_;
  return *this;
}

Vector3D &
Vector3D::operator-=(const Vector3D &v) {
  *impl_ -= *v.impl_;
  return *this;
}

Vector3D &
Vector3D::operator*=(double x) {
  *impl_ *= x;
  return *this;
}

Vector3D &
Vector3D::operator/=(double x) {
  *impl_ /= x;
  return *this;
}

double &
Vector3D::operator[](size_t i) {
  return (*impl_)[i];
}

const double &
Vector3D::operator[](size_t i) const {
  return (*impl_)[i];
}

Vector3D
Vector3D::operator-() const {
  return Vector3D(std::make_unique<Vector3DImpl>(-(*impl_)));
}

Vector3D
Vector3D::operator+(const Vector3D &v) const {
  return Vector3D(std::make_unique<Vector3DImpl>(*impl_ + *v.impl_));
}

Vector3D
Vector3D::operator-(const Vector3D &v) const {
  return Vector3D(std::make_unique<Vector3DImpl>(*impl_ - *v.impl_));
}

Vector3D
Vector3D::operator^(const Vector3D &v) const {
  return Vector3D(std::make_unique<Vector3DImpl>(*impl_ ^ *v.impl_));
}

double
Vector3D::operator*(const Vector3D &v) const {
  return *impl_ * *v.impl_;
}

Vector3D
Vector3D::operator*(double x) const {
  return Vector3D(std::make_unique<Vector3DImpl>(*impl_ * x));
}

Vector3D
Vector3D::operator/(double x) const {
  return Vector3D(std::make_unique<Vector3DImpl>(*impl_ / x));
}

double
Vector3D::norm() const {
  return impl_->norm();
}

double
Vector3D::normSqr() const {
  return impl_->normSqr();
}

Vector3D &
Vector3D::normalize() {
  impl_->normalize();
  return *this;
}

BSCurve::BSCurve()
  : impl_(std::make_unique<BSCurveImpl>()) {
}

BSCurve::BSCurve(size_t degree, DoubleVector knots, PointVector cpts)
  : impl_(std::make_unique<BSCurveImpl>(degree, knots, cpts)) {
}

BSCurve::BSCurve(const BSCurve &v)
  : impl_(std::make_unique<BSCurveImpl>(*v.impl_)) {
}

BSCurve::BSCurve(BSCurve &&v) = default;

BSCurve::~BSCurve() = default;

BSCurve &
BSCurve::operator=(BSCurve &&v) = default;

BSCurve &
BSCurve::operator=(const BSCurve &v) {
  *impl_ = *v.impl_;
  return *this;
}

Point3D
BSCurve::eval(double u) const {
  return impl_->eval(u);
}

Point3D
BSCurve::eval(double u, size_t nr_der, VectorVector &der) const {
  return impl_->eval(u, nr_der, der);
}

void
BSCurve::reverse() {
  impl_->reverse();
}

void
BSCurve::normalize() {
  impl_->normalize();
}

double
BSCurve::arcLength(double from, double to) const {
  return impl_->arcLength(from, to);
}

TriMesh::TriMesh()
  : impl_(std::make_unique<TriMeshImpl>()) {
}

TriMesh::TriMesh(const TriMesh &v)
  : impl_(std::make_unique<TriMeshImpl>(*v.impl_)) {
}

TriMesh::TriMesh(TriMesh &&v) = default;

TriMesh::~TriMesh() = default;

TriMesh &
TriMesh::operator=(TriMesh &&v) = default;

TriMesh &
TriMesh::operator=(const TriMesh &v) {
  *impl_ = *v.impl_;
  return *this;
}

void
TriMesh::resizePoints(size_t n) {
  impl_->resizePoints(n);
}

void
TriMesh::setPoint(size_t i, const Point3D &p) {
  impl_->setPoint(i, p);
}

void
TriMesh::setPoints(const PointVector &pv) {
  impl_->setPoints(pv);
}

void
TriMesh::addTriangle(size_t a, size_t b, size_t c) {
  impl_->addTriangle(a, b, c);
}

void
TriMesh::writeOBJ(std::string filename) const {
  impl_->writeOBJ(filename);
}
