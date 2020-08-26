#include <exception>

#include <Eigen/Sparse>

#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-dummy.hh"
#include "surface-harmonic.hh"
#include "utilities.hh"

namespace Transfinite {

using namespace Eigen;

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonDummy;

SurfaceHarmonic::SurfaceHarmonic() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceHarmonic::~SurfaceHarmonic() {
}

Point3D
SurfaceHarmonic::eval(const Point2D &uv) const {
  throw std::runtime_error("single-point evaluation is not supported for harmonic surfaces");
}

static Point3D
from2D(const Point2D &p) {
  return { p[0], p[1], 0 };
}

// Based on a program by Marton Vaitkus
TriMesh
SurfaceHarmonic::eval(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  size_t n_all = uvs.size(), n_boundary = 0;
  PointVector points; points.resize(n_all);

  // Set up valences
  std::vector<size_t> valences(n_all, 0);
  for (const auto &t : mesh.triangles())
    for (auto i : t)
      valences[i]++;

  // Fill the boundary points
  for (size_t i = 0; i < n_all; ++i)
    if (domain_->onEdge(resolution, i)) {
      auto bc = dynamic_cast<ParameterizationBarycentric*>(param_.get())->barycentric(uvs[i]);
      auto j = std::distance(bc.begin(), std::max_element(bc.begin(), bc.end()));
      if (bc[next(j)] > bc[prev(j)])
        j = next(j);
      points[i] = ribbons_[j]->curve()->eval(bc[j]);
      ++n_boundary;
    }

  // Set up the equations
  SparseMatrix<double> A(n_all + n_boundary, n_all + n_boundary);
  MatrixXd b = MatrixXd::Zero(n_all + n_boundary, 3);

  // - Sparse matrix handling
  VectorXd nnz(n_all + n_boundary);
  for (size_t i = 0; i < n_all; ++i)
    nnz(i) = valences[i] + 2;
  for (size_t i = 0; i < n_boundary; ++i)
    nnz(n_all + i) = 1;
  A.makeCompressed();
  A.reserve(nnz);

  // - Laplacian coeficients
  for (const auto &t : mesh.triangles()) {
    auto p1 = from2D(uvs[t[0]]), p2 = from2D(uvs[t[1]]), p3 = from2D(uvs[t[2]]);
    double Ai = ((p2 - p1) ^ (p3 - p1)).norm();
    double v1_cot = 0.5 * ((p2 - p1) * (p3 - p1)) / Ai;
    double v2_cot = 0.5 * ((p1 - p2) * (p3 - p2)) / Ai;
    double v3_cot = 0.5 * ((p2 - p3) * (p1 - p3)) / Ai;

    A.coeffRef(t[0], t[0]) += v3_cot + v2_cot;
    A.coeffRef(t[0], t[1]) += -v3_cot;
    A.coeffRef(t[0], t[2]) += -v2_cot;

    A.coeffRef(t[1], t[0]) += -v3_cot;
    A.coeffRef(t[1], t[1]) += v3_cot + v1_cot;
    A.coeffRef(t[1], t[2]) += -v1_cot;

    A.coeffRef(t[2], t[0]) += -v2_cot;
    A.coeffRef(t[2], t[1]) += -v1_cot;
    A.coeffRef(t[2], t[2]) += v2_cot + v1_cot;
  }

  // - Constraints
  for (size_t i = 0, j = 0; i < n_all; ++i)
    if (domain_->onEdge(resolution, i)) {
      A.coeffRef(n_all + j, i) = 1;
      A.coeffRef(i, n_all + j) = 1;
      b.block<1,3>(n_all + j, 0) = Map<const Vector3d>(points[i].data());
      ++j;
    }

  // Solve the system
  SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  MatrixXd x = solver.solve(b);

  for (size_t i = 0; i < n_all; ++i)
    points[i] = { x(i, 0), x(i, 1), x(i, 2) };

  mesh.setPoints(points);
  return mesh;
}

std::shared_ptr<Ribbon>
SurfaceHarmonic::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
