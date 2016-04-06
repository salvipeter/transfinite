#pragma once

#include <array>
#include <list>
#include <memory>
#include <vector>

namespace Geometry {

const double epsilon = 1.0e-8;

class Vector2D;
class Vector3D;
class BSCurve;

using Point2D = Vector2D;
using Point3D = Vector3D;
using DoubleVector = std::vector<double>;
using Vector2DVector = std::vector<Vector2D>;
using VectorVector = std::vector<Vector3D>;
using Point2DVector = std::vector<Point2D>;
using PointVector = std::vector<Point3D>;
using PointMatrix = std::vector<PointVector>;
using CurveVector = std::vector<std::shared_ptr<BSCurve>>;

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

class BSCurve {
public:
  // Constructors & destructor
  BSCurve();
  BSCurve(const PointVector &cpts);
  BSCurve(size_t degree, const DoubleVector &knots, const PointVector &cpts);
  BSCurve(const BSCurve &c);
  BSCurve(BSCurve &&c);
  ~BSCurve();

  // Assignments
  BSCurve &operator=(const BSCurve &c);
  BSCurve &operator=(BSCurve &&v);

  // Evaluation
  Point3D eval(double u) const;
  Point3D eval(double u, size_t nr_der, VectorVector &der) const;

  // Parameterization
  size_t degree() const;
  DoubleVector knotVector() const;
  void insertKnot(double k);
  void reverse();
  void normalize();

  // Control
  size_t nrControlPoints() const;
  Point3D controlPoint(size_t i) const;

  // Other
  double arcLength(double from, double to) const;
  void trim(double from, double to);

private:
  class BSCurveImpl;
  std::unique_ptr<BSCurveImpl> impl_;
};

struct BSSurface {
  size_t deg_u_, deg_v_;
  DoubleVector knots_u_, knots_v_;
  PointMatrix cnet_;
  CurveVector curves_, param_curves_;
};

class TriMesh {
public:
  using Triangle = std::array<size_t, 3>;

  // Constructors & destructor
  TriMesh();
  TriMesh(const TriMesh &c);
  TriMesh(TriMesh &&c);
  ~TriMesh();

  // Assignments
  TriMesh &operator=(const TriMesh &c);
  TriMesh &operator=(TriMesh &&v);

  // Mesh building
  void resizePoints(size_t n);
  void setPoints(const PointVector &pv);
  void addTriangle(size_t a, size_t b, size_t c);

  // I/O
  PointVector points() const;
  std::list<Triangle> triangles() const;
  Triangle closestTriangle(const Point3D &p) const;
  bool writeOBJ(std::string filename) const;

private:
  class TriMeshImpl;
  std::unique_ptr<TriMeshImpl> impl_;
};

class CurveFitter {
public:
  // Constructor & destructor
  CurveFitter();
  ~CurveFitter();

  // Fitting constraints
  void setDegree(size_t deg);
  void setNrControlPoints(size_t nr_cpts);
  void setKnotVector(const DoubleVector &knots);
  void addControlPoint(size_t i, const Point3D &point);
  void newPointGroup(double tolerance);
  void addParamPoint(double param, const Point3D &point);

  // Fit
  void fit();
  BSCurve curve() const;
  DoubleVector parameters(size_t group) const;

private:
  class CurveFitterImpl;
  std::unique_ptr<CurveFitterImpl> impl_;
};

class SurfaceFitter {
public:
  // Constructor & destructor
  SurfaceFitter();
  ~SurfaceFitter();

  // Fitting constraints
  void setDegreeU(size_t deg_u);
  void setDegreeV(size_t deg_v);
  void setNrControlPointsU(size_t nr_cpts_u);
  void setNrControlPointsV(size_t nr_cpts_v);
  void setKnotVectorU(const DoubleVector &knots_u);
  void setKnotVectorV(const DoubleVector &knots_v);
  void setMaxNrControlPointsU(size_t max_cpts_u);
  void setMaxNrControlPointsV(size_t max_cpts_v);
  void setCurvatureWeight(double weight);
  void setOscillationWeight(double weight);
  void addControlPoint(size_t i, size_t j, const Point3D &point);
  void newPointGroup(double tolerance);
  void addParamPoint(const Point2D &param, const Point3D &point);

  // Fit
  void fit();
  void fitWithCarrierSurface();
  BSSurface surface() const;
  Point2DVector parameters(size_t group) const;

private:
  class SurfaceFitterImpl;
  std::unique_ptr<SurfaceFitterImpl> impl_;
};

class IGES {
public:
  // Constructor & destructor
  IGES(std::string filename);
  ~IGES();

  // I/O
  void writeSurface(const BSSurface &surface);
  void close();

private:
  class IGESImpl;
  std::unique_ptr<IGESImpl> impl_;
};

} // namespace Geometry
