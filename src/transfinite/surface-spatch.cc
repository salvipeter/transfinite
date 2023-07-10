#include <algorithm>
#include <cmath>

#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-dummy.hh"
#include "surface-spatch.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonDummy;

SurfaceSPatch::SurfaceSPatch() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceSPatch::~SurfaceSPatch() {
}

static size_t multinomial(const SurfaceSPatch::Index &index) {
  auto fact = [](size_t n) { return (size_t)std::lround(std::tgamma(n + 1)); };
  size_t numerator = 0, denominator = 1;
  for (auto i : index) {
    numerator += i;
    denominator *= fact(i);
  }
  return fact(numerator) / denominator;
}

static double multiBernstein(const SurfaceSPatch::Index &index, const DoubleVector &bc) {
  double result = multinomial(index);
  for (size_t i = 0; i < index.size(); ++i)
    result *= std::pow(bc[i], index[i]);
  return result;
}

Point3D
SurfaceSPatch::eval(const Point2D &uv) const {
  DoubleVector bc =
    dynamic_cast<const ParameterizationBarycentric *>(param_.get())->barycentric(uv);
  Point3D p(0,0,0);
  for (const auto &cp : net_)
    p += cp.second * multiBernstein(cp.first, bc);
  return p;
}

void
SurfaceSPatch::initNetwork(size_t n, size_t d) {
  n_ = n;
  depth_ = d;
  net_.clear();
}

void
SurfaceSPatch::setupLoop() {
  CurveVector curves;
  for (size_t i = 0; i < n_; ++i) {
    size_t ip = (i + 1) % n_;
    PointVector pv;
    Index index(n_, 0);
    index[i] = depth_;
    for (size_t j = 0; j < depth_; ++j) {
      pv.push_back(net_[index]);
      --index[i];
      ++index[ip];
    }
    pv.push_back(net_[index]);
    curves.push_back(std::make_shared<BSplineCurve>(BSCurve(pv)));
  }
  setCurves(curves);

  Surface::setupLoop();

  if (domain_->update())
    param_->update();
}

void
SurfaceSPatch::setControlPoint(const Index &i, const Point3D &p) {
  net_[i] = p;
}

Point3D
SurfaceSPatch::controlPoint(const Index &i) const {
  return net_.at(i);
}

std::shared_ptr<Ribbon>
SurfaceSPatch::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
