#include <limits>

#include "vector-impl.hh"
#include "bspline-impl.hh"
#include "trimesh-impl.hh"
#include "fitter-impl.hh"
#include "iges-impl.hh"

namespace Geometry {

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

BSCurve::BSCurve(const PointVector &cpts)
  : impl_(std::make_unique<BSCurveImpl>(cpts)) {
}

BSCurve::BSCurve(size_t degree, const DoubleVector &knots, const PointVector &cpts)
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

size_t
BSCurve::degree() const {
  return impl_->degree();
}

DoubleVector
BSCurve::knotVector() const {
  return impl_->knotVector();
}

void
BSCurve::insertKnot(double k) {
  impl_->insertKnot(k);
}

void
BSCurve::reverse() {
  impl_->reverse();
}

void
BSCurve::normalize() {
  impl_->normalize();
}

size_t
BSCurve::nrControlPoints() const {
  return impl_->nrControlPoints();
}

Point3D
BSCurve::controlPoint(size_t i) const {
  return impl_->controlPoint(i);
}

double
BSCurve::arcLength(double from, double to) const {
  return impl_->arcLength(from, to);
}

void
BSCurve::trim(double from, double to) {
  impl_->trim(from, to);
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

PointVector
TriMesh::points() const {
  return impl_->points();
}

std::list<TriMesh::Triangle>
TriMesh::triangles() const {
  return impl_->triangles();
}

void
TriMesh::writeOBJ(std::string filename) const {
  impl_->writeOBJ(filename);
}

CurveFitter::CurveFitter()
  : impl_(std::make_unique<CurveFitterImpl>()) {
}

CurveFitter::~CurveFitter() = default;

void
CurveFitter::setDegree(size_t deg) {
  impl_->setDegree(deg);
}

void
CurveFitter::setNrControlPoints(size_t nr_cpts) {
  impl_->setNrControlPoints(nr_cpts);
}

void
CurveFitter::setKnotVector(const DoubleVector &knots) {
  impl_->setKnotVector(knots);
}

void
CurveFitter::addControlPoint(size_t i, const Point3D &point) {
  impl_->addControlPoint(i, point);
}

void
CurveFitter::newPointGroup(double tolerance) {
  impl_->newPointGroup(tolerance);
}

void
CurveFitter::addParamPoint(double param, const Point3D &point) {
  impl_->addParamPoint(param, point);
}

void
CurveFitter::fit() {
  impl_->fit();
}

BSCurve
CurveFitter::curve() const {
  return impl_->curve();
}

DoubleVector
CurveFitter::parameters(size_t group) const {
  return impl_->parameters(group);
}

SurfaceFitter::SurfaceFitter()
  : impl_(std::make_unique<SurfaceFitterImpl>()) {
}

SurfaceFitter::~SurfaceFitter() = default;

void
SurfaceFitter::setDegreeU(size_t deg_u) {
  impl_->setDegreeU(deg_u);
}

void
SurfaceFitter::setDegreeV(size_t deg_v) {
  impl_->setDegreeV(deg_v);
}

void
SurfaceFitter::setNrControlPointsU(size_t nr_cpts_u) {
  impl_->setNrControlPointsU(nr_cpts_u);
}

void
SurfaceFitter::setNrControlPointsV(size_t nr_cpts_v) {
  impl_->setNrControlPointsV(nr_cpts_v);
}

void
SurfaceFitter::setKnotVectorU(const DoubleVector &knots_u) {
  impl_->setKnotVectorU(knots_u);
}

void
SurfaceFitter::setKnotVectorV(const DoubleVector &knots_v) {
  impl_->setKnotVectorV(knots_v);
}

void
SurfaceFitter::setMaxNrControlPointsU(size_t max_cpts_u) {
  impl_->setMaxNrControlPointsU(max_cpts_u);
}

void
SurfaceFitter::setMaxNrControlPointsV(size_t max_cpts_v) {
  impl_->setMaxNrControlPointsV(max_cpts_v);
}

void
SurfaceFitter::setCurvatureWeight(double weight) {
  impl_->setCurvatureWeight(weight);
}

void
SurfaceFitter::setOscillationWeight(double weight) {
  impl_->setOscillationWeight(weight);
}

void
SurfaceFitter::addControlPoint(size_t i, size_t j, const Point3D &point) {
  impl_->addControlPoint(i, j, point);
}

void
SurfaceFitter::newPointGroup(double tolerance) {
  impl_->newPointGroup(tolerance);
}

void
SurfaceFitter::addParamPoint(const Point2D &param, const Point3D &point) {
  impl_->addParamPoint(param, point);
}

void
SurfaceFitter::fit() {
  impl_->fit();
}

void
SurfaceFitter::fitWithCarrierSurface() {
  impl_->fitWithCarrierSurface();
}

BSSurface
SurfaceFitter::surface() const {
  return impl_->surface();
}

Point2DVector
SurfaceFitter::parameters(size_t group) const {
  return impl_->parameters(group);
}

IGES::IGES(std::string filename)
  : impl_(std::make_unique<IGESImpl>(filename)) {
}

IGES::~IGES() = default;

void
IGES::writeSurface(const BSSurface &surface) {
  impl_->writeSurface(surface);
}

void
IGES::close() {
  impl_->close();
}

} // namespace Geometry
