#include <Eigen/Dense>

#include "domain-regular.hh"
#include "parameterization-bilinear.hh"
#include "ribbon-compatible-with-handler.hh"
#include "surface-elastic.hh"
#include "utilities.hh"

namespace Transfinite {

using namespace Eigen;

using DomainType = DomainRegular;
using ParamType = ParameterizationBilinear;
using RibbonType = RibbonCompatibleWithHandler;

SurfaceElastic::SurfaceElastic() : domain_tolerance(1.0e-5) {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceElastic::~SurfaceElastic() {
}

Point3D
SurfaceElastic::eval(const Point2D &uv) const {
  double u = uv[0], v = uv[1];
  auto sds = param_->mapToRibbons(uv);
  auto dpoly = domain()->vertices();

  Matrix3d A = Matrix3d::Zero(), b = Matrix3d::Zero();

  for (size_t i = 0; i < n_; ++i) {
    auto s = inrange(0, sds[i][0], 1);
    auto boundary_point = affineCombine(dpoly[prev(i)], s, dpoly[i]);
    auto ui = boundary_point[0], vi = boundary_point[1];
    auto di = (uv - boundary_point).norm();
    auto Pi = ribbons_[i]->eval({ s, 0 });
    if (di < domain_tolerance)
      return Pi;
    auto Ti = ribbons_[i]->crossDerivative(s) * di;

    double denom = std::pow(di, -3);
    A(0, 0) += 2 * denom;
    A(0, 1) += (u - ui) * denom;
    A(0, 2) += (v - vi) * denom;
    A(1, 0) += 3 * (u - ui) * denom;
    A(1, 1) += 2 * std::pow(u - ui, 2) * denom;
    A(1, 2) += 2 * (u - ui) * (v - vi) * denom;
    A(2, 0) += 3 * (v - vi) * denom;
    A(2, 1) += 2 * (u - ui) * (v - vi) * denom;
    A(2, 2) += 2 * std::pow(v - vi, 2) * denom;

    for (size_t k = 0; k < 3; ++k) {
      b(0, k) += (2 * Pi[k] + Ti[k]) * denom;
      b(1, k) += (3 * (u - ui) * Pi[k] + (u - ui) * Ti[k]) * denom;
      b(2, k) += (3 * (v - vi) * Pi[k] + (v - vi) * Ti[k]) * denom;
    }
  }

  Matrix3d x = A.colPivHouseholderQr().solve(b);
  Point3D p(x(0,0), x(0,1), x(0,2));
  // Vector3D j1(-x(1,0), -x(1,1), -x(1,2)), j2(-x(2,0), -x(2,1), -x(2,2));
  // auto normal_vector = (j1 ^ j2).normalize();
  return p;
}

std::shared_ptr<Ribbon>
SurfaceElastic::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
