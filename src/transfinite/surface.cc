#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"
#include "ribbon.hh"
#include "surface.hh"

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
Surface::setCurve(size_t i, const std::shared_ptr<BSCurve> &curve) {
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
  for (size_t i = 0; i < n_; ++i)
    ribbons_[i]->curve()->normalize();
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
      if (std::min(start_to_start, start_to_end) < std::min(end_to_start, end_to_end)) {
        ribbons_[i]->curve()->reverse();
        ribbons_[i]->curve()->normalize();
      }
    } else {
      Point3D r_start = ribbons_[i]->curve()->eval(0.0);
      Point3D r_end = ribbons_[i]->curve()->eval(1.0);
      Point3D rp_end = rp->curve()->eval(1.0);
      if ((r_end - rp_end).norm() < (r_start - rp_end).norm()) {
        ribbons_[i]->curve()->reverse();
        ribbons_[i]->curve()->normalize();
      }
    }
  }
}

void
Surface::update(size_t i) {
  if (domain_->update())
    param_->update();
  ribbons_[i]->update();
}

void
Surface::update() {
  if (domain_->update())
    param_->update();
  for (auto &r : ribbons_)
    r->update();
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

double
Surface::gamma(double d) const {
  if (use_gamma_)
    return d / (2.0 * d + 1.0);
  return d;
}
