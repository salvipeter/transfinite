#include <algorithm>
#include <exception>
#include <fstream>

#include "domain-regular.hh"
#include "parameterization-superd.hh"
#include "ribbon-dummy.hh"
#include "surface-superd.hh"
#include "utilities.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationSuperD;
using RibbonType = RibbonDummy;

SurfaceSuperD::SurfaceSuperD() : fullness_(0.5) {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceSuperD::~SurfaceSuperD() {
}

SurfaceSuperD::QuarticCurve
SurfaceSuperD::generateQuartic(const Point3D &a, const Point3D &b, const Point3D &c) const {
  double x1 = (2.0/5.0 * fullness_ + 3.0/5.0) * fullness_;
  double x2 = (-2.0/7.0 * fullness_ + 9.0/7.0) * fullness_;
  return
    { a,
      affineCombine(a, x1, b),
      affineCombine(affineCombine(a, x2, b), 0.5, affineCombine(c, x2, b)),
      affineCombine(c, x1, b),
      c };
}

SurfaceSuperD::QuarticCurve
SurfaceSuperD::generateBase(size_t i) const {
  return generateQuartic(cp_f_[prev(i)], cp_e_[i], cp_f_[i]);
}

SurfaceSuperD::QuarticCurve
SurfaceSuperD::generateMid(size_t i) const {
  return generateQuartic(cp_e_[prev(i)], cp_v_, cp_e_[next(i)]);
}

SurfaceSuperD::QuarticCurve
SurfaceSuperD::generateOpp(size_t i) const {
  if (n_ == 3) {
    QuarticCurve result;
    result.fill(cp_f_[next(i)]);
    return result;
  }
  if (n_ == 4) {
    auto result = generateBase(next(next(i)));
    std::reverse(result.begin(), result.end());
    return result;
  }
  throw std::runtime_error("generateOpp() should only be called for 3- and 4-sided patches");
}

// static
// void writeRibbonMesh(const SurfaceSuperD::QuarticSurface &s, const std::string &filename) {
//   std::ofstream f(filename);
//   for (size_t j = 0; j <= 4; ++j)
//     for (size_t k = 0; k <= 4; ++k) {
//       auto p = s[j][k];
//       f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
//     }
//   for (size_t j = 0; j <= 4; ++j)
//     for (size_t k = 0; k < 4; ++k) {
//       size_t base = j * 5 + k;
//       f << "l " << base + 1 << ' ' << base + 2 << std::endl;
//       base = k * 5 + j;
//       f << "l " << base + 1 << ' ' << base + 6 << std::endl;
//     }
// }

SurfaceSuperD::QuarticSurface
SurfaceSuperD::generateRibbon(size_t i) const {
  QuarticSurface result;
  auto base = generateBase(i);
  if (n_ > 4) {
    auto left = generateBase(prev(i));
    auto right = generateBase(next(i));
    for (size_t j = 0; j <= 4; ++j)
      for (size_t k = 0; k <= 4; ++k)
        result[j][k] =
          base[j] + (left[4-k] - left[4]) * (4 - j) / 4.0 + (right[k] - right[0]) * j / 4.0;
    return result;
  }
  auto mid = generateMid(i);
  auto opp = generateOpp(i);
  for (size_t j = 0; j <= 4; ++j)
    result[j] = generateQuartic(base[j], mid[j], opp[j]);
  return result;
}

static
Point3D bezierEvaluate(const SurfaceSuperD::QuarticSurface &s, double u, double v) {
  DoubleVector coeff_u, coeff_v;
  bernstein(4, u, coeff_u);
  bernstein(4, v, coeff_v);
  Point3D result(0, 0, 0);
  for (size_t i = 0; i <= 4; ++i)
    for (size_t j = 0; j <= 4; ++j)
      result += s[i][j] * coeff_u[i] * coeff_v[j];
  return result;
}

Point3D
SurfaceSuperD::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendSideSingular(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += bezierEvaluate(quartic_ribbons_[i], sds[i][0], sds[i][1]) * blends[i];
  return p;
}

void
SurfaceSuperD::initNetwork(size_t n) {
  n_ = n;
  cp_f_.resize(n);
  cp_e_.resize(n);
}

void
SurfaceSuperD::setupLoop() {
  CurveVector curves;
  for (size_t i = 0; i < n_; ++i) {
    auto base = generateBase(i);
    PointVector pv;
    pv.assign(base.begin(), base.end());
    curves.push_back(std::make_shared<BSplineCurve>(BSCurve(pv)));
  }
  setCurves(curves);

  Surface::setupLoop();

  if (domain_->update())
    param_->update();
}

void
SurfaceSuperD::updateRibbons() {
  quartic_ribbons_.clear();
  for (size_t i = 0; i < n_; ++i)
    quartic_ribbons_.emplace_back(generateRibbon(i));
}

Point3D
SurfaceSuperD::vertexControlPoint() const {
  return cp_v_;
}

void
SurfaceSuperD::setVertexControlPoint(const Point3D &p) {
  cp_v_ = p;
}

Point3D
SurfaceSuperD::faceControlPoint(size_t i) const {
  return cp_f_[i];
}

void
SurfaceSuperD::setFaceControlPoint(size_t i, const Point3D &p) {
  cp_f_[i] = p;
}

Point3D
SurfaceSuperD::edgeControlPoint(size_t i) const {
  return cp_e_[i];
}

void
SurfaceSuperD::setEdgeControlPoint(size_t i, const Point3D &p) {
  cp_e_[i] = p;
}

double
SurfaceSuperD::fullness() const {
  return fullness_;
}

void
SurfaceSuperD::setFullness(double f) {
  fullness_ = f;
}

std::shared_ptr<Ribbon>
SurfaceSuperD::newRibbon() const {
  return std::make_shared<RibbonType>();
}

} // namespace Transfinite
