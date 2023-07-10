#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"
#include "ribbon.hh"
#include "surface.hh"
#include "utilities.hh"

// For the 4-sided central split fit
#include "parameterization-bilinear.hh"

namespace Transfinite {

Surface::Surface()
  : n_(0), use_gamma_(true) {
}

Surface::~Surface() {
}

void
Surface::setGamma(bool use_gamma) {
  use_gamma_ = use_gamma;
}

void
Surface::setCurve(size_t i, const std::shared_ptr<Curve> &curve) {
  if (n_ <= i) {
    ribbons_.resize(i + 1);
    n_ = i + 1;
  }
  ribbons_[i] = newRibbon();
  ribbons_[i]->setCurve(curve);
  domain_->setSide(i, curve);
}

void
Surface::setCurves(const CurveVector &curves) {
  ribbons_.clear();
  ribbons_.reserve(curves.size());
  for (const auto &c : curves) {
    ribbons_.push_back(newRibbon());
    ribbons_.back()->setCurve(c);
  }
  domain_->setSides(curves);
  n_ = curves.size();
}

void
Surface::setupLoop() {
  // Tasks:
  // - propagate adjacency information
  // - normalize curves
  // - reverse curves when needed (and normalize once again, for safety)
  for (size_t i = 0; i < n_; ++i) {
    std::shared_ptr<Ribbon> rp = ribbons_[prev(i)], rn = ribbons_[next(i)];
    ribbons_[i]->setNeighbors(rp, rn);
    if (i == 0) {
      Point3D r_start = ribbons_[i]->curve()->eval(0.0);
      Point3D r_end = ribbons_[i]->curve()->eval(1.0);
      Point3D rn_start = rn->curve()->eval(0.0);
      Point3D rn_end = rn->curve()->eval(1.0);
      double end_to_start = (r_end - rn_start).norm();
      double end_to_end = (r_end - rn_end).norm();
      double start_to_start = (r_start - rn_start).norm();
      double start_to_end = (r_start - rn_end).norm();
      if (std::min(start_to_start, start_to_end) < std::min(end_to_start, end_to_end))
        ribbons_[i]->curve()->reverse();
    } else {
      Point3D r_start = ribbons_[i]->curve()->eval(0.0);
      Point3D r_end = ribbons_[i]->curve()->eval(1.0);
      Point3D rp_end = rp->curve()->eval(1.0);
      if ((r_end - rp_end).norm() < (r_start - rp_end).norm())
        ribbons_[i]->curve()->reverse();
    }
  }
}

std::optional<Vector3D>
Surface::ribbonHandler(size_t i) const {
  return ribbons_[i]->handler();
}

void
Surface::setRibbonHandler(size_t i, const Vector3D &h) {
  ribbons_[i]->setHandler(h);
}

void
Surface::overrideNormalFence(size_t i, const std::shared_ptr<NormalFence> &fence) {
  ribbons_[i]->overrideNormalFence(fence);
}

double
Surface::ribbonMultiplier(size_t i) const {
  return ribbons_[i]->multiplier();
}

void
Surface::setRibbonMultiplier(size_t i, double m) {
  ribbons_[i]->setMultiplier(m);
}

void
Surface::resetRibbon(size_t i) {
  ribbons_[i]->reset();
}

void
Surface::update(size_t i) {
  if (domain_->update())
    param_->update();
  ribbons_[i]->update();
  updateCorner(prev(i));
  updateCorner(i);
}

void
Surface::update() {
  if (domain_->update())
    param_->update();
  for (auto &r : ribbons_)
    r->update();
  updateCorners();
}

std::shared_ptr<const Domain>
Surface::domain() const {
  return domain_;
}

std::shared_ptr<const Parameterization>
Surface::parameterization() const {
  return param_;
}

std::shared_ptr<const Ribbon>
Surface::ribbon(size_t i) const {
  return ribbons_[i];
}

TriMesh
Surface::eval(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [this](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}

Point3D
Surface::cornerCorrection(size_t i, double s1, double s2) const {
  // Assumes that both s1 and s2 are 0 at the corner,
  // s1 increases towards corner (i-1), and s2 towards corner (i+1).
  s1 = inrange(0, gamma(s1), 1);
  s2 = inrange(0, gamma(s2), 1);
  return corner_data_[i].point
    + corner_data_[i].tangent1 * s1
    + corner_data_[i].tangent2 * s2
    + rationalTwist(s1, s2, corner_data_[i].twist2, corner_data_[i].twist1) * s1 * s2;
}

Point3D
Surface::sideInterpolant(size_t i, double si, double di) const {
  si = inrange(0, si, 1);
  di = std::max(gamma(di), 0.0);
  return ribbons_[i]->eval(Point2D(si, di));
}

Point3D
Surface::cornerInterpolant(size_t i, const Point2DVector &sds) const {
  double si = sds[i][0], si1 = sds[next(i)][0];
  return sideInterpolant(i, si, si1) + sideInterpolant(next(i), si1, 1.0 - si)
    - cornerCorrection(i, 1.0 - si, si1);
}

Point3D
Surface::cornerInterpolantD(size_t i, const Point2DVector &sds) const {
  double di = sds[i][1], di1 = sds[next(i)][1];
  return sideInterpolant(i, 1.0 - di1, di) + sideInterpolant(next(i), di, di1)
    - cornerCorrection(i, di1, di);
}

DoubleVector
Surface::blendCorner(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);

  size_t close_to_boundary = 0;
  for (const auto &sd : sds) {
    if (sd[1] < epsilon)
      ++close_to_boundary;
  }

  if (close_to_boundary > 0) {
    for (size_t i = 0; i < n_; ++i) {
      size_t ip = next(i);
      if (close_to_boundary > 1)
        blf.push_back(sds[i][1] < epsilon && sds[ip][1] < epsilon ? 1.0 : 0.0);
      else if (sds[i][1] < epsilon) {
        double tmp = std::pow(sds[ip][1], -2);
        blf.push_back(tmp / (tmp + std::pow(sds[prev(i)][1], -2)));
      } else if (sds[ip][1] < epsilon) {
        double tmp = std::pow(sds[i][1], -2);
        blf.push_back(tmp / (tmp + std::pow(sds[next(ip)][1], -2)));
      } else
        blf.push_back(0.0);
    }
  } else {
    double denominator = 0.0;
    for (size_t i = 0; i < n_; ++i) {
      blf.push_back(std::pow(sds[i][1] * sds[next(i)][1], -2));
      denominator += blf.back();
    }
    std::transform(blf.begin(), blf.end(), blf.begin(),
                   [denominator](double x) { return x / denominator; });
  }

  return blf;
}

DoubleVector
Surface::blendSideSingular(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);

  size_t close_to_boundary = 0;
  for (const auto &sd : sds) {
    if (sd[1] < epsilon)
      ++close_to_boundary;
  }

  if (close_to_boundary > 0) {
    double blend_val = 1.0 / close_to_boundary;
    for (const auto &sd : sds)
      blf.push_back(sd[1] < epsilon ? blend_val : 0.0);
  } else {
    double denominator = 0.0;
    for (const auto &sd : sds) {
      blf.push_back(std::pow(sd[1], -2));
      denominator += blf.back();
    }
    std::transform(blf.begin(), blf.end(), blf.begin(),
                   [denominator](double x) { return x / denominator; });
  }

  return blf;
}

DoubleVector
Surface::blendCornerDeficient(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    size_t ip = next(i);
    if (sds[i][1] < epsilon && sds[ip][1] < epsilon) {
      blf.push_back(1.0);
      continue;
    }
    blf.push_back((sds[ip][1] * hermite(0, 1.0 - sds[i][0]) * hermite(0, sds[i][1] ) +
                   sds[i][1]  * hermite(0,    sds[ip][0]  ) * hermite(0, sds[ip][1])) /
                  (sds[i][1] + sds[ip][1]));
  }
  return blf;
}

void
Surface::updateCorner(size_t i) {
  static const double step = 1.0e-4;
  size_t ip = next(i);

  VectorVector der;
  Vector3D d1, d2;
  corner_data_[i].point = ribbons_[i]->curve()->eval(1.0, 1, der);
  corner_data_[i].tangent1 = -der[1];
  ribbons_[ip]->curve()->eval(0.0, 1, der);
  corner_data_[i].tangent2 = der[1];
  d1 = ribbons_[i]->crossDerivative(1.0);
  d2 = ribbons_[i]->crossDerivative(1.0 - step);
  corner_data_[i].twist1 = (d2 - d1) / step;
  d1 = ribbons_[ip]->crossDerivative(0.0);
  d2 = ribbons_[ip]->crossDerivative(step);
  corner_data_[i].twist2 = (d2 - d1) / step;
}

void
Surface::updateCorners() {
  corner_data_.resize(n_);
  for (size_t i = 0; i < n_; ++i)
    updateCorner(i);
}

double
Surface::gamma(double d) const {
  if (use_gamma_)
    return d / (2.0 * d + 1.0);
  return d;
}

Vector3D
Surface::rationalTwist(double u, double v, const Vector3D &f, const Vector3D &g) {
  if (std::abs(u + v) < epsilon)
    return Vector3D(0,0,0);
  return (f * u + g * v) / (u + v);
}

} // namespace Transfinite
