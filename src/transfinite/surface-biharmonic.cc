#include <exception>

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include "Eigen/SVD"

#ifdef HAVE_LIBTRIANGLE
#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle.h>
}
#endif // HAVE_LIBTRIANGLE

#include "domain-angular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-compatible.hh"
#include "surface-biharmonic.hh"
#include "utilities.hh"

// Based on:
// T. Stanko, S. Hahmann, G-P Bonneau, N. Saguin-Sprynski:
//   Surfacing Curve Networks with Normal Control.
//     Computers and Graphics, Vol. 60, pp. 1-8, 2016.

namespace Transfinite {

using namespace Eigen;

using DomainType = DomainAngular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonCompatible;

SurfaceBiharmonic::SurfaceBiharmonic() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceBiharmonic::~SurfaceBiharmonic() {
}

Point3D
SurfaceBiharmonic::eval(const Point2D &uv) const {
  throw std::runtime_error("single-point evaluation is not supported for biharmonic surfaces");
}

static SparseMatrix<double> laplaceMatrix(const TriMesh &mesh, const PointVector &uvs) {
  size_t n_all = uvs.size();

  // Set up valences
  std::vector<size_t> valences(n_all, 0);
  for (const auto &t : mesh.triangles())
    for (auto i : t)
      valences[i]++;

  // Sparse matrix handling
  SparseMatrix<double> Ls(n_all, n_all);
  VectorXd nnz(n_all);
  for (size_t i = 0; i < n_all; ++i)
    nnz(i) = valences[i] + 1;
  Ls.makeCompressed();
  Ls.reserve(nnz);

  // Laplacian coeficients
  for (const auto &t : mesh.triangles()) {
    auto p1 = uvs[t[0]], p2 = uvs[t[1]], p3 = uvs[t[2]];
    auto a = p3 - p2, b = p1 - p3, c = p2 - p1;
    double Ai = (b ^ c).norm();
    double v1_cot = -0.5 * (c * b) / Ai;
    double v2_cot = -0.5 * (c * a) / Ai;
    double v3_cot = -0.5 * (a * b) / Ai;

    Ls.coeffRef(t[0], t[0]) += v3_cot + v2_cot;
    Ls.coeffRef(t[0], t[1]) += -v3_cot;
    Ls.coeffRef(t[0], t[2]) += -v2_cot;

    Ls.coeffRef(t[1], t[0]) += -v3_cot;
    Ls.coeffRef(t[1], t[1]) += v3_cot + v1_cot;
    Ls.coeffRef(t[1], t[2]) += -v1_cot;

    Ls.coeffRef(t[2], t[0]) += -v2_cot;
    Ls.coeffRef(t[2], t[1]) += -v1_cot;
    Ls.coeffRef(t[2], t[2]) += v2_cot + v1_cot;
  }

  return Ls;
}

static double computeAngle(Vector3D u, Vector3D v) {
  u.normalize(); v.normalize();
  return std::acos(std::min(std::max(u * v, -1.0), 1.0));
}

static double voronoiArea(const Point3D &p1, const Point3D &p2, const Point3D &p3) {
  double a2 = (p3 - p2).normSqr(), b2 = (p1 - p3).normSqr(), c2 = (p2 - p1).normSqr();
  double alpha = computeAngle(p2 - p1, p3 - p1);

  if (a2 + b2 < c2)                // obtuse gamma
    return 0.125 * b2 * std::tan(alpha);
  if (a2 + c2 < b2)                // obtuse beta
    return 0.125 * c2 * std::tan(alpha);
  if (b2 + c2 < a2) {              // obtuse alpha
    double b = std::sqrt(b2), c = std::sqrt(c2);
    double total_area = 0.5 * b * c * std::sin(alpha);
    double beta = computeAngle(p1 - p2, p3 - p2);
    double gamma = computeAngle(p1 - p3, p2 - p3);
    return total_area - 0.125 * (b2 * std::tan(gamma) + c2 * std::tan(beta));
  }

  double r2 = 0.25 * a2 / std::pow(std::sin(alpha), 2); // squared circumradius
  auto area = [r2](double x2) {
    return 0.125 * std::sqrt(x2) * std::sqrt(std::max(4.0 * r2 - x2, 0.0));
  };
  return area(b2) + area(c2);
}

static DoubleVector voronoiAreas(const TriMesh &mesh, const PointVector &uvs) {
  size_t n_all = uvs.size();
  DoubleVector areas(n_all);
  for (const auto &t : mesh.triangles()) {
    auto p1 = uvs[t[0]], p2 = uvs[t[1]], p3 = uvs[t[2]];
    areas[t[0]] += voronoiArea(p1, p2, p3);
    areas[t[1]] += voronoiArea(p2, p3, p1);
    areas[t[2]] += voronoiArea(p3, p1, p2);
  }  
  return areas;
}

static DoubleVector computeCurvatures(const TriMesh &mesh, const DoubleVector &areas,
                                      std::function<bool(size_t i)> on_edge,
                                      const MatrixXd &points, const MatrixXd &normals) {
  size_t n_all = areas.size();
  MatrixXd mean = MatrixXd::Zero(n_all, 3);
  for (const auto &t : mesh.triangles()) {
    auto addComponent = [&](size_t i1, size_t i2, size_t i3) {
      const Vector3d &p1 = points.row(i1), &p2 = points.row(i2), &p3 = points.row(i3);
      const Vector3d &n1 = normals.row(i1), &n2 = normals.row(i2), &n3 = normals.row(i3);
      if (on_edge(i1)) {
        if (on_edge(i2))
          mean.row(i1) += (n1 * 2 + n2).normalized().cross(p2 - p1);
        if (on_edge(i3))
          mean.row(i1) += (n1 * 2 + n3).normalized().cross(p1 - p3);
      }
      mean.row(i1) += (n1 + n2 + n3).normalized().cross(p3 - p2);
    };
    addComponent(t[0], t[1], t[2]);
    addComponent(t[1], t[2], t[0]);
    addComponent(t[2], t[0], t[1]);
  }
  DoubleVector result(n_all);
  for (size_t i = 0; i < n_all; ++i)
    result[i] = mean.row(i).norm() / (2 * areas[i]);
  return result;
}

static SparseMatrix<double> prepareSolver(const TriMesh &mesh, const PointVector &points,
                                          const std::vector<size_t> &boundary, bool propagation) {
  size_t n_all = points.size(), n_boundary = boundary.size();

  // Compute the required matrices
  SparseMatrix<double> Ls = laplaceMatrix(mesh, points);
  auto areas = voronoiAreas(mesh, points);
  Map<VectorXd> M(&areas[0], areas.size());
  DiagonalMatrix<double,Dynamic> M1 = M.asDiagonal().inverse();
  for (size_t b : boundary)
    M1.diagonal()(b) = 0;               // hack (private correspondence with T. Stanko)
  SparseMatrix<double> L = M1 * Ls;

  // Set up the equations
  SparseMatrix<double> A = propagation ? Ls * L : SparseMatrix<double>(L.transpose() * L);
  A.conservativeResize(n_all + n_boundary, n_all + n_boundary);

  // Add the lagrange multiplicator part
  for (size_t j = 0; j < n_boundary; ++j) {
    size_t i = boundary[j];
    A.coeffRef(n_all + j, i) = 1;
    A.coeffRef(i, n_all + j) = 1;
  }

  return A;
}

[[maybe_unused]]
static Point2DVector projectToLSQPlane(const PointVector &pv) {
  // Compute centroid
  size_t n = pv.size();
  Point3D centroid(0, 0, 0);
  for (const auto &p : pv)
    centroid += p;
  centroid /= n;

  // Solve by singular value decomposition
  Eigen::MatrixXd A(n, 3);
  for (size_t i = 0; i < n; ++i) {
    auto p = pv[i] - centroid;
    A.row(i) << p[0], p[1], p[2];
  }
  auto svd = A.jacobiSvd(Eigen::ComputeThinV);

  // Extract principal directions
  auto vec = [](const auto &v) { return Vector3D(v(0), v(1), v(2)); };
  auto u = vec(svd.matrixV().col(0));
  auto v = vec(svd.matrixV().col(1));

  Point2DVector result;
  for (const auto &p : pv)
    result.emplace_back((p - centroid) * u, (p - centroid) * v);
  return result;
}

[[maybe_unused]]
auto
SurfaceBiharmonic::generateDomainOld(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  auto on_edge = [=](size_t i) { return domain_->onEdge(resolution, i); };
  auto vertex_boundary = [=](size_t i) {
    auto bc = dynamic_cast<ParameterizationBarycentric*>(param_.get())->barycentric(uvs[i]);
    auto j = std::distance(bc.begin(), std::max_element(bc.begin(), bc.end()));
    if (bc[next(j)] > bc[prev(j)])
      j = next(j);
    return std::make_pair(j, bc[j]);
  };
  return std::make_tuple(mesh, uvs, on_edge, vertex_boundary);
}

[[maybe_unused]]
auto
SurfaceBiharmonic::generateDomain(size_t resolution) const {
#ifdef HAVE_LIBTRIANGLE
  PointVector points3d;
  double length = 0;
  for (size_t i = 0; i < domain_->size(); ++i) {
    auto curve = ribbons_[i]->curve();
    length += curve->arcLength(curve->basis().low(), curve->basis().high());
    for (size_t j = 0; j < resolution; ++j) {
      double u = (double)j / resolution;
      points3d.push_back(curve->eval(u));
    }
  }
  length /= resolution * domain_->size();

  auto projected = projectToLSQPlane(points3d);
  double max_area = (length * length * std::sqrt(3.0)) / 4;

  // Input points
  size_t n = projected.size();
  DoubleVector points; points.reserve(2 * n);
  for (const auto &p : projected) {
    points.push_back(p[0]);
    points.push_back(p[1]);
  }

  // Input segments : just a closed polygon
  std::vector<int> segments; segments.reserve(2 * n);
  for (size_t i = 0; i < n; ++i) {
    segments.push_back(i);
    segments.push_back(i + 1);
  }
  segments.back() = 0;
    
  // Setup output data structure
  struct triangulateio in, out;
  in.pointlist = &points[0];
  in.numberofpoints = n;
  in.numberofpointattributes = 0;
  in.pointmarkerlist = nullptr;
  in.segmentlist = &segments[0];
  in.numberofsegments = n;
  in.segmentmarkerlist = nullptr;
  in.numberofholes = 0;
  in.numberofregions = 0;

  // Setup output data structure
  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.triangleattributelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  // Call the library function [with maximum triangle area = resolution]
  std::ostringstream cmd;
  cmd << "pqa" << std::fixed << max_area << "DBPzQY";
  triangulate(const_cast<char *>(cmd.str().c_str()), &in, &out, (struct triangulateio *)nullptr);

  TriMesh mesh;
  Point2DVector uvs;
  for (int i = 0; i < out.numberofpoints; ++i)
    uvs.emplace_back(out.pointlist[2*i], out.pointlist[2*i+1]);
  for (int i = 0; i < out.numberoftriangles; ++i)
    mesh.addTriangle(out.trianglelist[3*i+0],
                     out.trianglelist[3*i+1],
                     out.trianglelist[3*i+2]);

  auto on_edge = [=](size_t i) { return i < n; };
  auto vertex_boundary = [=](size_t i) {
    size_t j = i / resolution;
    double u = (double)(i - j * resolution) / resolution;
    return std::make_pair(j, u);
  };
  return std::make_tuple(mesh, uvs, on_edge, vertex_boundary);
#else  // !HAVE_LIBTRIANGLE
  return generateDomainOld(resolution);
#endif // HAVE_LIBTRIANGLE
}

TriMesh
SurfaceBiharmonic::eval(size_t resolution) const {
  auto [mesh, uvs, on_edge, vertex_boundary] = generateDomain(resolution);
  PointVector uvs3d;
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(uvs3d),
                 [](const Point2D &p) { return Point3D(p[0], p[1], 0); });

  size_t n_all = uvs.size();
  PointVector points(n_all);
  VectorVector normals(n_all);

  // Set up boundary index map
  std::vector<size_t> boundary;
  for (size_t i = 0; i < n_all; ++i)
    if (on_edge(i))
      boundary.push_back(i);
  size_t n_boundary = boundary.size();

  // Fill the boundary points & normals
  for (auto i : boundary) {
    auto [j, u] = vertex_boundary(i);
    points[i] = ribbons_[j]->curve()->eval(u);
    normals[i] = ribbons_[j]->normal(u);
  }

  // Compute the system
  auto A = prepareSolver(mesh, uvs3d, boundary, true);
  SparseLU<SparseMatrix<double>> solver(A);

  // Propagate the boundary points
  MatrixXd b = MatrixXd::Zero(n_all + n_boundary, 3);
  for (size_t j = 0; j < n_boundary; ++j)
    b.row(n_all + j) = Map<const Vector3d>(points[boundary[j]].data());
  MatrixXd Vstar = solver.solve(b);

  // Propagate the normal vectors
  b = MatrixXd::Zero(n_all + n_boundary, 3);
  for (size_t j = 0; j < n_boundary; ++j)
    b.row(n_all + j) = Map<const Vector3d>(normals[boundary[j]].data());
  MatrixXd Nstar = solver.solve(b);
  for (size_t i = 0; i < n_all; ++i)
    Nstar.row(i) = Nstar.row(i).normalized();

  // Recompute the matrices with the propagated surface
  for (size_t i = 0; i < n_all; ++i)
    points[i] = { Vstar(i, 0), Vstar(i, 1), Vstar(i, 2) };
  auto A2 = prepareSolver(mesh, points, boundary, false);
  SparseLU<SparseMatrix<double>> solver2(A2);

  // Also recompute L & areas (TODO: redundant)
  SparseMatrix<double> Ls = laplaceMatrix(mesh, points);
  auto areas = voronoiAreas(mesh, points);
  Map<VectorXd> M(&areas[0], areas.size());
  SparseMatrix<double> L = M.asDiagonal().inverse() * Ls;

  // Compute mean curvature
  auto mean = computeCurvatures(mesh, areas, on_edge, Vstar, Nstar);
  MatrixXd H(n_all, 3);
  for (size_t i = 0; i < n_all; ++i)
    H.row(i) = - Nstar.row(i) * mean[i];
  MatrixXd LH = L * H;

  // Compute the final surface
  b = MatrixXd::Zero(n_all + n_boundary, 3);
  b.block(0, 0, n_all, 3) = LH;
  for (size_t j = 0; j < n_boundary; ++j)
    b.row(n_all + j) = Vstar.row(boundary[j]);
  MatrixXd x = solver2.solve(b);

  for (size_t i = 0; i < n_all; ++i)
    points[i] = { x(i, 0), x(i, 1), x(i, 2) };

  mesh.setPoints(points);
  return mesh;
}

std::shared_ptr<Ribbon>
SurfaceBiharmonic::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
